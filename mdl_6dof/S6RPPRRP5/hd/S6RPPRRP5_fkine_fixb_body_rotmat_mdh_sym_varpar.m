% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRP5 (for one body)
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
% Datum: 2020-11-04 21:28
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:28:55
	% EndTime: 2020-11-04 21:28:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:28:55
	% EndTime: 2020-11-04 21:28:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t45 = cos(qJ(1));
	t44 = sin(qJ(1));
	t1 = [t45, -t44, 0, 0; t44, t45, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:28:55
	% EndTime: 2020-11-04 21:28:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t1 = [0, -t47, t46, t47 * pkin(1) + t46 * qJ(2) + 0; 0, -t46, -t47, t46 * pkin(1) - t47 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:28:55
	% EndTime: 2020-11-04 21:28:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->4)
	t50 = cos(qJ(1));
	t49 = sin(qJ(1));
	t48 = pkin(1) + qJ(3);
	t1 = [0, t49, t50, t49 * qJ(2) + t48 * t50 + 0; 0, -t50, t49, -t50 * qJ(2) + t48 * t49 + 0; 1, 0, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:28:55
	% EndTime: 2020-11-04 21:28:55
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (13->11), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->7)
	t56 = cos(qJ(1));
	t55 = cos(qJ(4));
	t54 = sin(qJ(1));
	t53 = sin(qJ(4));
	t52 = pkin(1) + qJ(3);
	t51 = pkin(7) - qJ(2);
	t1 = [t56 * t53, t56 * t55, -t54, -t51 * t54 + t52 * t56 + 0; t54 * t53, t54 * t55, t56, t51 * t56 + t52 * t54 + 0; t55, -t53, 0, pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:28:55
	% EndTime: 2020-11-04 21:28:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (24->20), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t59 = sin(qJ(5));
	t61 = sin(qJ(1));
	t68 = t61 * t59;
	t62 = cos(qJ(5));
	t67 = t61 * t62;
	t64 = cos(qJ(1));
	t66 = t64 * t59;
	t65 = t64 * t62;
	t63 = cos(qJ(4));
	t60 = sin(qJ(4));
	t58 = pkin(7) - qJ(2);
	t57 = t60 * pkin(4) - t63 * pkin(8) + pkin(1) + qJ(3);
	t1 = [t60 * t65 - t68, -t60 * t66 - t67, -t64 * t63, t57 * t64 - t58 * t61 + 0; t60 * t67 + t66, -t60 * t68 + t65, -t61 * t63, t57 * t61 + t58 * t64 + 0; t63 * t62, -t63 * t59, t60, t63 * pkin(4) + t60 * pkin(8) + pkin(2) + pkin(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:28:55
	% EndTime: 2020-11-04 21:28:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (32->23), mult. (31->22), div. (0->0), fcn. (44->6), ass. (0->15)
	t73 = sin(qJ(5));
	t75 = sin(qJ(1));
	t82 = t75 * t73;
	t76 = cos(qJ(5));
	t81 = t75 * t76;
	t78 = cos(qJ(1));
	t80 = t78 * t73;
	t79 = t78 * t76;
	t77 = cos(qJ(4));
	t74 = sin(qJ(4));
	t72 = qJ(6) + pkin(8);
	t71 = t76 * pkin(5) + pkin(4);
	t70 = t73 * pkin(5) + pkin(7) - qJ(2);
	t69 = t71 * t74 - t72 * t77 + pkin(1) + qJ(3);
	t1 = [t74 * t79 - t82, -t74 * t80 - t81, -t78 * t77, t69 * t78 - t70 * t75 + 0; t74 * t81 + t80, -t74 * t82 + t79, -t75 * t77, t69 * t75 + t70 * t78 + 0; t77 * t76, -t77 * t73, t74, t77 * t71 + t72 * t74 + pkin(2) + pkin(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end