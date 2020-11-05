% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR8 (for one body)
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

function Tc_mdh = S5RPPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:35
	% EndTime: 2020-11-04 20:15:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:35
	% EndTime: 2020-11-04 20:15:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t43 = cos(qJ(1));
	t42 = sin(qJ(1));
	t1 = [t43, -t42, 0, 0; t42, t43, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:35
	% EndTime: 2020-11-04 20:15:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t45 = cos(qJ(1));
	t44 = sin(qJ(1));
	t1 = [t45, 0, t44, t45 * pkin(1) + t44 * qJ(2) + 0; t44, 0, -t45, t44 * pkin(1) - t45 * qJ(2) + 0; 0, 1, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:35
	% EndTime: 2020-11-04 20:15:35
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t52 = pkin(1) + pkin(2);
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t49 = cos(pkin(8));
	t48 = sin(pkin(8));
	t47 = -t51 * t48 + t50 * t49;
	t46 = -t50 * t48 - t51 * t49;
	t1 = [-t46, t47, 0, t50 * qJ(2) + t52 * t51 + 0; t47, t46, 0, -t51 * qJ(2) + t52 * t50 + 0; 0, 0, -1, -qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:35
	% EndTime: 2020-11-04 20:15:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (25->14), mult. (16->10), div. (0->0), fcn. (24->6), ass. (0->10)
	t61 = cos(qJ(1));
	t60 = sin(qJ(1));
	t59 = pkin(8) + qJ(4);
	t58 = cos(t59);
	t57 = sin(t59);
	t56 = sin(pkin(8)) * pkin(3) + qJ(2);
	t55 = cos(pkin(8)) * pkin(3) + pkin(1) + pkin(2);
	t54 = -t61 * t57 + t60 * t58;
	t53 = -t60 * t57 - t61 * t58;
	t1 = [-t53, t54, 0, t55 * t61 + t56 * t60 + 0; t54, t53, 0, t55 * t60 - t56 * t61 + 0; 0, 0, -1, -pkin(6) - qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:35
	% EndTime: 2020-11-04 20:15:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (46->23), mult. (26->18), div. (0->0), fcn. (40->10), ass. (0->15)
	t69 = pkin(8) + qJ(4);
	t75 = pkin(1) + pkin(2);
	t74 = cos(qJ(1));
	t73 = cos(qJ(5));
	t72 = sin(qJ(1));
	t71 = sin(qJ(5));
	t70 = -qJ(1) + pkin(8);
	t68 = -qJ(1) + t69;
	t67 = cos(t69);
	t66 = sin(t69);
	t65 = cos(t68);
	t64 = sin(t68);
	t63 = t72 * t66 + t74 * t67;
	t62 = t74 * t66 - t72 * t67;
	t1 = [t63 * t73, -t63 * t71, t62, pkin(4) * t65 + pkin(7) * t64 + pkin(3) * cos(t70) + t75 * t74 + t72 * qJ(2) + 0; -t62 * t73, t62 * t71, t63, pkin(7) * t65 - pkin(4) * t64 - pkin(3) * sin(t70) + t75 * t72 - t74 * qJ(2) + 0; -t71, -t73, 0, -pkin(6) - qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end