% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:15
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:51
	% EndTime: 2020-11-04 20:15:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:51
	% EndTime: 2020-11-04 20:15:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t42 = cos(qJ(1));
	t41 = sin(qJ(1));
	t1 = [t42, -t41, 0, 0; t41, t42, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:51
	% EndTime: 2020-11-04 20:15:51
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t44 = cos(qJ(1));
	t43 = sin(qJ(1));
	t1 = [t44, 0, t43, t44 * pkin(1) + t43 * qJ(2) + 0; t43, 0, -t44, t43 * pkin(1) - t44 * qJ(2) + 0; 0, 1, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:51
	% EndTime: 2020-11-04 20:15:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t51 = pkin(1) + pkin(2);
	t50 = cos(qJ(1));
	t49 = sin(qJ(1));
	t48 = cos(pkin(8));
	t47 = sin(pkin(8));
	t46 = -t50 * t47 + t49 * t48;
	t45 = -t49 * t47 - t50 * t48;
	t1 = [-t45, t46, 0, t49 * qJ(2) + t51 * t50 + 0; t46, t45, 0, -t50 * qJ(2) + t51 * t49 + 0; 0, 0, -1, -qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:51
	% EndTime: 2020-11-04 20:15:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (26->17), mult. (28->16), div. (0->0), fcn. (42->6), ass. (0->11)
	t62 = cos(qJ(1));
	t61 = cos(qJ(4));
	t60 = sin(qJ(1));
	t59 = sin(qJ(4));
	t58 = cos(pkin(8));
	t57 = sin(pkin(8));
	t55 = -t57 * pkin(3) + t58 * pkin(6) - qJ(2);
	t54 = t58 * pkin(3) + t57 * pkin(6) + pkin(1) + pkin(2);
	t53 = t60 * t57 + t62 * t58;
	t52 = t62 * t57 - t60 * t58;
	t1 = [t53 * t61, -t53 * t59, t52, t54 * t62 - t55 * t60 + 0; -t52 * t61, t52 * t59, t53, t54 * t60 + t55 * t62 + 0; -t59, -t61, 0, -qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:51
	% EndTime: 2020-11-04 20:15:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (41->23), mult. (66->30), div. (0->0), fcn. (89->8), ass. (0->16)
	t78 = sin(qJ(1));
	t70 = sin(qJ(5));
	t73 = cos(qJ(4));
	t77 = t70 * t73;
	t72 = cos(qJ(5));
	t76 = t72 * t73;
	t71 = sin(qJ(4));
	t75 = pkin(4) * t73 + pkin(7) * t71 + pkin(3);
	t74 = cos(qJ(1));
	t69 = cos(pkin(8));
	t68 = sin(pkin(8));
	t66 = t78 * t68 + t74 * t69;
	t65 = t74 * t68 - t78 * t69;
	t64 = -t69 * pkin(6) + t75 * t68 + qJ(2);
	t63 = t68 * pkin(6) + t75 * t69 + pkin(1) + pkin(2);
	t1 = [t65 * t70 + t66 * t76, t72 * t65 - t66 * t77, t66 * t71, t63 * t74 + t64 * t78 + 0; -t65 * t76 + t66 * t70, t65 * t77 + t66 * t72, -t65 * t71, t63 * t78 - t64 * t74 + 0; -t71 * t72, t71 * t70, t73, -t71 * pkin(4) + t73 * pkin(7) + pkin(5) - qJ(3) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end