% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPPR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:11
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:49
	% EndTime: 2020-11-04 20:11:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:49
	% EndTime: 2020-11-04 20:11:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t50 = cos(qJ(1));
	t49 = sin(qJ(1));
	t1 = [t50, -t49, 0, 0; t49, t50, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:49
	% EndTime: 2020-11-04 20:11:49
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t54 = cos(qJ(1));
	t53 = sin(qJ(1));
	t52 = cos(pkin(7));
	t51 = sin(pkin(7));
	t1 = [t54 * t52, -t54 * t51, t53, pkin(1) * t54 + qJ(2) * t53 + 0; t53 * t52, -t53 * t51, -t54, pkin(1) * t53 - qJ(2) * t54 + 0; t51, t52, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:49
	% EndTime: 2020-11-04 20:11:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->14), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t59 = cos(qJ(1));
	t58 = sin(qJ(1));
	t57 = cos(pkin(7));
	t56 = sin(pkin(7));
	t55 = pkin(2) * t57 + t56 * qJ(3) + pkin(1);
	t1 = [t58, -t59 * t57, t59 * t56, t58 * qJ(2) + t55 * t59 + 0; -t59, -t58 * t57, t58 * t56, -t59 * qJ(2) + t55 * t58 + 0; 0, -t56, -t57, t56 * pkin(2) - t57 * qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:49
	% EndTime: 2020-11-04 20:11:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->14)
	t61 = sin(pkin(8));
	t67 = sin(qJ(1));
	t72 = t67 * t61;
	t63 = cos(pkin(8));
	t71 = t67 * t63;
	t68 = cos(qJ(1));
	t70 = t68 * t61;
	t69 = t68 * t63;
	t66 = pkin(2) + qJ(4);
	t65 = pkin(3) + qJ(2);
	t64 = cos(pkin(7));
	t62 = sin(pkin(7));
	t60 = t62 * qJ(3) + t66 * t64 + pkin(1);
	t1 = [t62 * t70 + t71, t62 * t69 - t72, t68 * t64, t60 * t68 + t65 * t67 + 0; t62 * t72 - t69, t62 * t71 + t70, t67 * t64, t60 * t67 - t65 * t68 + 0; -t64 * t61, -t64 * t63, t62, -t64 * qJ(3) + t66 * t62 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:49
	% EndTime: 2020-11-04 20:11:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (38->26), mult. (57->36), div. (0->0), fcn. (78->8), ass. (0->21)
	t78 = sin(pkin(7));
	t84 = cos(qJ(5));
	t92 = t78 * t84;
	t80 = cos(pkin(7));
	t91 = t80 * t84;
	t82 = sin(qJ(5));
	t90 = t82 * t80;
	t77 = sin(pkin(8));
	t83 = sin(qJ(1));
	t89 = t83 * t77;
	t79 = cos(pkin(8));
	t88 = t83 * t79;
	t85 = cos(qJ(1));
	t87 = t85 * t77;
	t86 = t85 * t79;
	t81 = pkin(2) + qJ(4);
	t76 = -t77 * pkin(4) + t79 * pkin(6) - qJ(3);
	t75 = t79 * pkin(4) + t77 * pkin(6) + pkin(3) + qJ(2);
	t74 = t78 * t87 + t88;
	t73 = -t76 * t78 + t81 * t80 + pkin(1);
	t1 = [t74 * t84 + t85 * t90, -t74 * t82 + t85 * t91, -t78 * t86 + t89, t73 * t85 + t75 * t83 + 0; (t77 * t92 + t90) * t83 - t84 * t86, (-t78 * t89 + t86) * t82 + t83 * t91, -t78 * t88 - t87, t73 * t83 - t75 * t85 + 0; -t77 * t91 + t78 * t82, t77 * t90 + t92, t80 * t79, t76 * t80 + t81 * t78 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end