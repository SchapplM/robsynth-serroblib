% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRP6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:29
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:15
	% EndTime: 2020-11-04 21:29:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:15
	% EndTime: 2020-11-04 21:29:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t1 = [t47, -t46, 0, 0; t46, t47, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:15
	% EndTime: 2020-11-04 21:29:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t49 = cos(qJ(1));
	t48 = sin(qJ(1));
	t1 = [0, -t49, t48, t49 * pkin(1) + t48 * qJ(2) + 0; 0, -t48, -t49, t48 * pkin(1) - t49 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:15
	% EndTime: 2020-11-04 21:29:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->4)
	t52 = cos(qJ(1));
	t51 = sin(qJ(1));
	t50 = pkin(1) + qJ(3);
	t1 = [0, t51, t52, t51 * qJ(2) + t50 * t52 + 0; 0, -t52, t51, -t52 * qJ(2) + t50 * t51 + 0; 1, 0, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:15
	% EndTime: 2020-11-04 21:29:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (13->11), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->7)
	t58 = cos(qJ(1));
	t57 = cos(qJ(4));
	t56 = sin(qJ(1));
	t55 = sin(qJ(4));
	t54 = pkin(1) + qJ(3);
	t53 = pkin(7) - qJ(2);
	t1 = [t58 * t55, t58 * t57, -t56, -t53 * t56 + t54 * t58 + 0; t56 * t55, t56 * t57, t58, t53 * t58 + t54 * t56 + 0; t57, -t55, 0, pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:15
	% EndTime: 2020-11-04 21:29:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (24->20), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t61 = sin(qJ(5));
	t63 = sin(qJ(1));
	t70 = t63 * t61;
	t64 = cos(qJ(5));
	t69 = t63 * t64;
	t66 = cos(qJ(1));
	t68 = t66 * t61;
	t67 = t66 * t64;
	t65 = cos(qJ(4));
	t62 = sin(qJ(4));
	t60 = pkin(7) - qJ(2);
	t59 = t62 * pkin(4) - t65 * pkin(8) + pkin(1) + qJ(3);
	t1 = [t62 * t67 - t70, -t62 * t68 - t69, -t66 * t65, t59 * t66 - t60 * t63 + 0; t62 * t69 + t68, -t62 * t70 + t67, -t63 * t65, t59 * t63 + t60 * t66 + 0; t65 * t64, -t65 * t61, t62, t65 * pkin(4) + t62 * pkin(8) + pkin(2) + pkin(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:29:15
	% EndTime: 2020-11-04 21:29:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (33->23), mult. (36->24), div. (0->0), fcn. (49->6), ass. (0->14)
	t74 = sin(qJ(5));
	t76 = sin(qJ(1));
	t83 = t76 * t74;
	t77 = cos(qJ(5));
	t82 = t76 * t77;
	t79 = cos(qJ(1));
	t81 = t79 * t74;
	t80 = t79 * t77;
	t78 = cos(qJ(4));
	t75 = sin(qJ(4));
	t73 = pkin(5) * t77 + qJ(6) * t74 + pkin(4);
	t72 = -t74 * pkin(5) + qJ(6) * t77 - pkin(7) + qJ(2);
	t71 = -t78 * pkin(8) + t73 * t75 + pkin(1) + qJ(3);
	t1 = [t75 * t80 - t83, -t79 * t78, t75 * t81 + t82, t71 * t79 + t72 * t76 + 0; t75 * t82 + t81, -t76 * t78, t75 * t83 - t80, t71 * t76 - t72 * t79 + 0; t78 * t77, t75, t78 * t74, t75 * pkin(8) + t73 * t78 + pkin(2) + pkin(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end