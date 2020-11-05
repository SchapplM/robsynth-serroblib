% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:27
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:27:15
	% EndTime: 2020-11-04 20:27:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:27:15
	% EndTime: 2020-11-04 20:27:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t1 = [t47, -t46, 0, 0; t46, t47, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:27:15
	% EndTime: 2020-11-04 20:27:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t50 = qJ(1) + pkin(9);
	t49 = cos(t50);
	t48 = sin(t50);
	t1 = [t49, -t48, 0, cos(qJ(1)) * pkin(1) + 0; t48, t49, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:27:15
	% EndTime: 2020-11-04 20:27:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t55 = cos(qJ(3));
	t54 = sin(qJ(3));
	t53 = qJ(1) + pkin(9);
	t52 = cos(t53);
	t51 = sin(t53);
	t1 = [t52 * t55, -t52 * t54, t51, t52 * pkin(2) + t51 * pkin(6) + cos(qJ(1)) * pkin(1) + 0; t51 * t55, -t51 * t54, -t52, t51 * pkin(2) - t52 * pkin(6) + sin(qJ(1)) * pkin(1) + 0; t54, t55, 0, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:27:15
	% EndTime: 2020-11-04 20:27:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (38->19), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->11)
	t59 = sin(qJ(4));
	t62 = cos(qJ(3));
	t65 = t59 * t62;
	t61 = cos(qJ(4));
	t64 = t61 * t62;
	t60 = sin(qJ(3));
	t63 = pkin(3) * t62 + pkin(7) * t60 + pkin(2);
	t58 = qJ(1) + pkin(9);
	t57 = cos(t58);
	t56 = sin(t58);
	t1 = [t56 * t59 + t57 * t64, t56 * t61 - t57 * t65, t57 * t60, cos(qJ(1)) * pkin(1) + t56 * pkin(6) + 0 + t63 * t57; t56 * t64 - t57 * t59, -t56 * t65 - t57 * t61, t56 * t60, sin(qJ(1)) * pkin(1) - t57 * pkin(6) + 0 + t63 * t56; t60 * t61, -t60 * t59, -t62, t60 * pkin(3) - t62 * pkin(7) + pkin(5) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:27:15
	% EndTime: 2020-11-04 20:27:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (58->23), mult. (39->26), div. (0->0), fcn. (52->10), ass. (0->15)
	t72 = qJ(4) + qJ(5);
	t69 = sin(t72);
	t75 = cos(qJ(3));
	t80 = t69 * t75;
	t70 = cos(t72);
	t79 = t70 * t75;
	t78 = pkin(4) * sin(qJ(4)) + pkin(6);
	t66 = cos(qJ(4)) * pkin(4) + pkin(3);
	t74 = sin(qJ(3));
	t76 = -pkin(8) - pkin(7);
	t77 = t66 * t75 - t74 * t76 + pkin(2);
	t71 = qJ(1) + pkin(9);
	t68 = cos(t71);
	t67 = sin(t71);
	t1 = [t67 * t69 + t68 * t79, t67 * t70 - t68 * t80, t68 * t74, cos(qJ(1)) * pkin(1) + 0 + t78 * t67 + t77 * t68; t67 * t79 - t68 * t69, -t67 * t80 - t68 * t70, t67 * t74, sin(qJ(1)) * pkin(1) + 0 - t78 * t68 + t77 * t67; t74 * t70, -t74 * t69, -t75, t74 * t66 + t75 * t76 + pkin(5) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end