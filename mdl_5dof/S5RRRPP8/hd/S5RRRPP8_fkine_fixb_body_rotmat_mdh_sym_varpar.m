% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPP8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:41
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRPP8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:41:35
	% EndTime: 2020-11-04 20:41:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:41:35
	% EndTime: 2020-11-04 20:41:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t43 = cos(qJ(1));
	t42 = sin(qJ(1));
	t1 = [t43, -t42, 0, 0; t42, t43, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:41:35
	% EndTime: 2020-11-04 20:41:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t47 = cos(qJ(1));
	t46 = cos(qJ(2));
	t45 = sin(qJ(1));
	t44 = sin(qJ(2));
	t1 = [t47 * t46, -t47 * t44, t45, t47 * pkin(1) + t45 * pkin(6) + 0; t45 * t46, -t45 * t44, -t47, t45 * pkin(1) - t47 * pkin(6) + 0; t44, t46, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:41:35
	% EndTime: 2020-11-04 20:41:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t51 = sin(qJ(1));
	t53 = cos(qJ(2));
	t57 = t51 * t53;
	t49 = sin(qJ(3));
	t54 = cos(qJ(1));
	t56 = t54 * t49;
	t52 = cos(qJ(3));
	t55 = t54 * t52;
	t50 = sin(qJ(2));
	t48 = t53 * pkin(2) + t50 * pkin(7) + pkin(1);
	t1 = [t51 * t49 + t53 * t55, t51 * t52 - t53 * t56, t54 * t50, t51 * pkin(6) + t48 * t54 + 0; t52 * t57 - t56, -t49 * t57 - t55, t51 * t50, -t54 * pkin(6) + t48 * t51 + 0; t50 * t52, -t50 * t49, -t53, t50 * pkin(2) - t53 * pkin(7) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:41:35
	% EndTime: 2020-11-04 20:41:35
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (27->19), mult. (36->25), div. (0->0), fcn. (49->6), ass. (0->13)
	t63 = sin(qJ(1));
	t65 = cos(qJ(2));
	t69 = t63 * t65;
	t61 = sin(qJ(3));
	t66 = cos(qJ(1));
	t68 = t66 * t61;
	t64 = cos(qJ(3));
	t67 = t66 * t64;
	t62 = sin(qJ(2));
	t60 = -t61 * pkin(3) + qJ(4) * t64 - pkin(6);
	t59 = t64 * pkin(3) + t61 * qJ(4) + pkin(2);
	t58 = t62 * pkin(7) + t59 * t65 + pkin(1);
	t1 = [t66 * t62, -t63 * t61 - t65 * t67, -t63 * t64 + t65 * t68, t58 * t66 - t60 * t63 + 0; t63 * t62, -t64 * t69 + t68, t61 * t69 + t67, t58 * t63 + t60 * t66 + 0; -t65, -t62 * t64, t62 * t61, -t65 * pkin(7) + t59 * t62 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:41:35
	% EndTime: 2020-11-04 20:41:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (34->20), mult. (36->25), div. (0->0), fcn. (49->6), ass. (0->15)
	t75 = sin(qJ(1));
	t77 = cos(qJ(2));
	t83 = t75 * t77;
	t73 = sin(qJ(3));
	t78 = cos(qJ(1));
	t82 = t78 * t73;
	t76 = cos(qJ(3));
	t81 = t78 * t76;
	t72 = qJ(5) + pkin(3);
	t80 = qJ(4) * t76 - t72 * t73 - pkin(6);
	t79 = pkin(4) + pkin(7);
	t74 = sin(qJ(2));
	t71 = t73 * qJ(4) + t72 * t76 + pkin(2);
	t70 = t71 * t77 + t79 * t74 + pkin(1);
	t1 = [t78 * t74, -t75 * t76 + t77 * t82, t75 * t73 + t77 * t81, t70 * t78 - t80 * t75 + 0; t75 * t74, t73 * t83 + t81, t76 * t83 - t82, t70 * t75 + t80 * t78 + 0; -t77, t74 * t73, t74 * t76, t71 * t74 - t79 * t77 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end