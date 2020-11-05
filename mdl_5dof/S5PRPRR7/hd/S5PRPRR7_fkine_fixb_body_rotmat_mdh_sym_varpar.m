% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:01
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:01:07
	% EndTime: 2020-11-04 20:01:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:01:07
	% EndTime: 2020-11-04 20:01:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t49 = cos(pkin(8));
	t48 = sin(pkin(8));
	t1 = [t49, -t48, 0, 0; t48, t49, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:01:07
	% EndTime: 2020-11-04 20:01:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t53 = cos(qJ(2));
	t52 = sin(qJ(2));
	t51 = cos(pkin(8));
	t50 = sin(pkin(8));
	t1 = [t51 * t53, -t51 * t52, t50, t51 * pkin(1) + t50 * pkin(5) + 0; t50 * t53, -t50 * t52, -t51, t50 * pkin(1) - t51 * pkin(5) + 0; t52, t53, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:01:07
	% EndTime: 2020-11-04 20:01:07
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (16->14), mult. (18->12), div. (0->0), fcn. (26->4), ass. (0->6)
	t56 = sin(qJ(2));
	t57 = cos(qJ(2));
	t58 = pkin(2) * t57 + qJ(3) * t56 + pkin(1);
	t55 = cos(pkin(8));
	t54 = sin(pkin(8));
	t1 = [t54, -t55 * t57, t55 * t56, t54 * pkin(5) + t58 * t55 + 0; -t55, -t54 * t57, t54 * t56, -t55 * pkin(5) + t58 * t54 + 0; 0, -t56, -t57, t56 * pkin(2) - t57 * qJ(3) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:01:07
	% EndTime: 2020-11-04 20:01:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->17), mult. (30->22), div. (0->0), fcn. (43->6), ass. (0->12)
	t61 = sin(qJ(4));
	t62 = sin(qJ(2));
	t69 = t61 * t62;
	t63 = cos(qJ(4));
	t68 = t62 * t63;
	t64 = cos(qJ(2));
	t66 = pkin(2) + pkin(6);
	t67 = qJ(3) * t62 + t64 * t66 + pkin(1);
	t65 = pkin(3) + pkin(5);
	t60 = cos(pkin(8));
	t59 = sin(pkin(8));
	t1 = [t59 * t63 + t60 * t69, -t59 * t61 + t60 * t68, t60 * t64, t65 * t59 + t67 * t60 + 0; t59 * t69 - t60 * t63, t59 * t68 + t60 * t61, t59 * t64, t67 * t59 - t65 * t60 + 0; -t64 * t61, -t64 * t63, t62, -t64 * qJ(3) + t66 * t62 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:01:08
	% EndTime: 2020-11-04 20:01:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (40->20), mult. (37->24), div. (0->0), fcn. (50->8), ass. (0->14)
	t75 = sin(pkin(8));
	t77 = sin(qJ(2));
	t84 = t75 * t77;
	t76 = cos(pkin(8));
	t83 = t76 * t77;
	t82 = pkin(4) * cos(qJ(4)) + pkin(3) + pkin(5);
	t70 = sin(qJ(4)) * pkin(4) + qJ(3);
	t73 = pkin(2) + pkin(6) + pkin(7);
	t79 = cos(qJ(2));
	t81 = t70 * t77 + t73 * t79 + pkin(1);
	t74 = qJ(4) + qJ(5);
	t72 = cos(t74);
	t71 = sin(t74);
	t1 = [t71 * t83 + t75 * t72, -t75 * t71 + t72 * t83, t76 * t79, t82 * t75 + t81 * t76 + 0; t71 * t84 - t76 * t72, t76 * t71 + t72 * t84, t75 * t79, t81 * t75 - t82 * t76 + 0; -t79 * t71, -t79 * t72, t77, -t70 * t79 + t73 * t77 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end