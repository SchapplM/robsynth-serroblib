% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:32
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:32:02
	% EndTime: 2020-11-04 21:32:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:32:02
	% EndTime: 2020-11-04 21:32:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t44 = cos(qJ(1));
	t43 = sin(qJ(1));
	t1 = [t44, -t43, 0, 0; t43, t44, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:32:02
	% EndTime: 2020-11-04 21:32:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t46 = cos(qJ(1));
	t45 = sin(qJ(1));
	t1 = [0, -t46, t45, t46 * pkin(1) + t45 * qJ(2) + 0; 0, -t45, -t46, t45 * pkin(1) - t46 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:32:02
	% EndTime: 2020-11-04 21:32:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->4)
	t49 = cos(qJ(1));
	t48 = sin(qJ(1));
	t47 = pkin(1) + qJ(3);
	t1 = [0, t48, t49, t48 * qJ(2) + t47 * t49 + 0; 0, -t49, t48, -t49 * qJ(2) + t47 * t48 + 0; 1, 0, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:32:02
	% EndTime: 2020-11-04 21:32:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->7)
	t55 = cos(qJ(1));
	t54 = cos(qJ(4));
	t53 = sin(qJ(1));
	t52 = sin(qJ(4));
	t51 = pkin(1) + qJ(3);
	t50 = pkin(7) - qJ(2);
	t1 = [t55 * t52, t55 * t54, -t53, -t50 * t53 + t51 * t55 + 0; t53 * t52, t53 * t54, t55, t50 * t55 + t51 * t53 + 0; t54, -t52, 0, pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:32:02
	% EndTime: 2020-11-04 21:32:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (24->20), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t58 = sin(qJ(5));
	t60 = sin(qJ(1));
	t67 = t60 * t58;
	t61 = cos(qJ(5));
	t66 = t60 * t61;
	t63 = cos(qJ(1));
	t65 = t63 * t58;
	t64 = t63 * t61;
	t62 = cos(qJ(4));
	t59 = sin(qJ(4));
	t57 = pkin(7) - qJ(2);
	t56 = t59 * pkin(4) - t62 * pkin(8) + pkin(1) + qJ(3);
	t1 = [t59 * t64 - t67, -t59 * t65 - t66, -t63 * t62, t56 * t63 - t57 * t60 + 0; t59 * t66 + t65, -t59 * t67 + t64, -t60 * t62, t56 * t60 + t57 * t63 + 0; t62 * t61, -t62 * t58, t59, t62 * pkin(4) + t59 * pkin(8) + pkin(2) + pkin(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:32:02
	% EndTime: 2020-11-04 21:32:02
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (42->24), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t73 = qJ(5) + qJ(6);
	t71 = sin(t73);
	t75 = sin(qJ(1));
	t82 = t75 * t71;
	t72 = cos(t73);
	t81 = t75 * t72;
	t77 = cos(qJ(1));
	t80 = t77 * t71;
	t79 = t77 * t72;
	t78 = pkin(9) + pkin(8);
	t76 = cos(qJ(4));
	t74 = sin(qJ(4));
	t70 = cos(qJ(5)) * pkin(5) + pkin(4);
	t69 = sin(qJ(5)) * pkin(5) + pkin(7) - qJ(2);
	t68 = t70 * t74 - t76 * t78 + pkin(1) + qJ(3);
	t1 = [t74 * t79 - t82, -t74 * t80 - t81, -t77 * t76, t68 * t77 - t69 * t75 + 0; t74 * t81 + t80, -t74 * t82 + t79, -t75 * t76, t68 * t75 + t69 * t77 + 0; t76 * t72, -t76 * t71, t74, t70 * t76 + t74 * t78 + pkin(2) + pkin(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end