% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR11 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:44
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRPR11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:12
	% EndTime: 2020-11-04 20:44:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:12
	% EndTime: 2020-11-04 20:44:12
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t49 = cos(qJ(1));
	t48 = sin(qJ(1));
	t1 = [t49, -t48, 0, 0; t48, t49, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:12
	% EndTime: 2020-11-04 20:44:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t53 = cos(qJ(1));
	t52 = cos(qJ(2));
	t51 = sin(qJ(1));
	t50 = sin(qJ(2));
	t1 = [t53 * t52, -t53 * t50, t51, t53 * pkin(1) + t51 * pkin(6) + 0; t51 * t52, -t51 * t50, -t53, t51 * pkin(1) - t53 * pkin(6) + 0; t50, t52, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:12
	% EndTime: 2020-11-04 20:44:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t57 = sin(qJ(1));
	t59 = cos(qJ(2));
	t63 = t57 * t59;
	t55 = sin(qJ(3));
	t60 = cos(qJ(1));
	t62 = t60 * t55;
	t58 = cos(qJ(3));
	t61 = t60 * t58;
	t56 = sin(qJ(2));
	t54 = t59 * pkin(2) + t56 * pkin(7) + pkin(1);
	t1 = [t57 * t55 + t59 * t61, t57 * t58 - t59 * t62, t60 * t56, t57 * pkin(6) + t54 * t60 + 0; t58 * t63 - t62, -t55 * t63 - t61, t57 * t56, -t60 * pkin(6) + t54 * t57 + 0; t56 * t58, -t56 * t55, -t59, t56 * pkin(2) - t59 * pkin(7) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:13
	% EndTime: 2020-11-04 20:44:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (26->18), mult. (36->25), div. (0->0), fcn. (49->6), ass. (0->13)
	t69 = sin(qJ(1));
	t71 = cos(qJ(2));
	t75 = t69 * t71;
	t67 = sin(qJ(3));
	t72 = cos(qJ(1));
	t74 = t72 * t67;
	t70 = cos(qJ(3));
	t73 = t72 * t70;
	t68 = sin(qJ(2));
	t66 = -t67 * pkin(3) + qJ(4) * t70 - pkin(6);
	t65 = t70 * pkin(3) + t67 * qJ(4) + pkin(2);
	t64 = t68 * pkin(7) + t65 * t71 + pkin(1);
	t1 = [t69 * t67 + t71 * t73, t72 * t68, -t69 * t70 + t71 * t74, t64 * t72 - t66 * t69 + 0; t70 * t75 - t74, t69 * t68, t67 * t75 + t73, t64 * t69 + t66 * t72 + 0; t68 * t70, -t71, t68 * t67, -t71 * pkin(7) + t65 * t68 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:44:13
	% EndTime: 2020-11-04 20:44:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (45->23), mult. (56->30), div. (0->0), fcn. (79->8), ass. (0->18)
	t83 = sin(qJ(1));
	t86 = cos(qJ(2));
	t92 = t83 * t86;
	t87 = cos(qJ(1));
	t91 = t86 * t87;
	t81 = sin(qJ(3));
	t85 = cos(qJ(3));
	t89 = pkin(3) + pkin(4);
	t90 = qJ(4) * t85 - t89 * t81 - pkin(6);
	t88 = pkin(7) - pkin(8);
	t84 = cos(qJ(5));
	t82 = sin(qJ(2));
	t80 = sin(qJ(5));
	t79 = t81 * qJ(4) + t89 * t85 + pkin(2);
	t78 = t80 * t81 + t84 * t85;
	t77 = -t80 * t85 + t84 * t81;
	t76 = t79 * t86 + t88 * t82 + pkin(1);
	t1 = [t83 * t77 + t78 * t91, t77 * t91 - t83 * t78, -t87 * t82, t76 * t87 - t90 * t83 + 0; -t77 * t87 + t78 * t92, t77 * t92 + t78 * t87, -t83 * t82, t76 * t83 + t90 * t87 + 0; t82 * t78, t82 * t77, t86, t79 * t82 - t88 * t86 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end