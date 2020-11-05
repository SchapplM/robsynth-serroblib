% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:26
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:26:04
	% EndTime: 2020-11-04 21:26:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:26:04
	% EndTime: 2020-11-04 21:26:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t42 = cos(qJ(1));
	t41 = sin(qJ(1));
	t1 = [t42, -t41, 0, 0; t41, t42, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:26:04
	% EndTime: 2020-11-04 21:26:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t44 = cos(qJ(1));
	t43 = sin(qJ(1));
	t1 = [0, -t44, t43, t44 * pkin(1) + t43 * qJ(2) + 0; 0, -t43, -t44, t43 * pkin(1) - t44 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:26:04
	% EndTime: 2020-11-04 21:26:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->4)
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t45 = pkin(1) + qJ(3);
	t1 = [0, t46, t47, t46 * qJ(2) + t45 * t47 + 0; 0, -t47, t46, -t47 * qJ(2) + t45 * t46 + 0; 1, 0, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:26:04
	% EndTime: 2020-11-04 21:26:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->7)
	t53 = cos(qJ(1));
	t52 = cos(qJ(4));
	t51 = sin(qJ(1));
	t50 = sin(qJ(4));
	t49 = pkin(1) + qJ(3);
	t48 = pkin(7) - qJ(2);
	t1 = [t53 * t50, t53 * t52, -t51, -t48 * t51 + t49 * t53 + 0; t51 * t50, t51 * t52, t53, t48 * t53 + t49 * t51 + 0; t52, -t50, 0, pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:26:04
	% EndTime: 2020-11-04 21:26:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (24->20), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t55 = sin(pkin(9));
	t59 = sin(qJ(1));
	t65 = t59 * t55;
	t56 = cos(pkin(9));
	t64 = t59 * t56;
	t61 = cos(qJ(1));
	t63 = t61 * t55;
	t62 = t61 * t56;
	t60 = cos(qJ(4));
	t58 = sin(qJ(4));
	t57 = pkin(7) - qJ(2);
	t54 = t58 * pkin(4) - t60 * qJ(5) + pkin(1) + qJ(3);
	t1 = [t58 * t62 - t65, -t58 * t63 - t64, -t61 * t60, t54 * t61 - t57 * t59 + 0; t58 * t64 + t63, -t58 * t65 + t62, -t59 * t60, t54 * t59 + t57 * t61 + 0; t60 * t56, -t60 * t55, t58, t60 * pkin(4) + t58 * qJ(5) + pkin(2) + pkin(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:26:04
	% EndTime: 2020-11-04 21:26:04
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (42->24), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t71 = pkin(9) + qJ(6);
	t69 = sin(t71);
	t74 = sin(qJ(1));
	t80 = t74 * t69;
	t70 = cos(t71);
	t79 = t74 * t70;
	t76 = cos(qJ(1));
	t78 = t76 * t69;
	t77 = t76 * t70;
	t75 = cos(qJ(4));
	t73 = sin(qJ(4));
	t72 = qJ(5) + pkin(8);
	t68 = cos(pkin(9)) * pkin(5) + pkin(4);
	t67 = sin(pkin(9)) * pkin(5) + pkin(7) - qJ(2);
	t66 = t68 * t73 - t72 * t75 + pkin(1) + qJ(3);
	t1 = [t73 * t77 - t80, -t73 * t78 - t79, -t76 * t75, t66 * t76 - t67 * t74 + 0; t73 * t79 + t78, -t73 * t80 + t77, -t74 * t75, t66 * t74 + t67 * t76 + 0; t75 * t70, -t75 * t69, t73, t68 * t75 + t72 * t73 + pkin(2) + pkin(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end