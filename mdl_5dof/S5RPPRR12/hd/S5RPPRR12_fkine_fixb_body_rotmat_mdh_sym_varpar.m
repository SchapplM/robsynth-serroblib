% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR12 (for one body)
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
% Datum: 2020-11-04 20:16
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:40
	% EndTime: 2020-11-04 20:16:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:40
	% EndTime: 2020-11-04 20:16:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t42 = cos(qJ(1));
	t41 = sin(qJ(1));
	t1 = [t42, -t41, 0, 0; t41, t42, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:40
	% EndTime: 2020-11-04 20:16:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t44 = cos(qJ(1));
	t43 = sin(qJ(1));
	t1 = [0, -t44, t43, t44 * pkin(1) + t43 * qJ(2) + 0; 0, -t43, -t44, t43 * pkin(1) - t44 * qJ(2) + 0; 1, 0, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:40
	% EndTime: 2020-11-04 20:16:40
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t49 = cos(qJ(1));
	t48 = sin(qJ(1));
	t47 = pkin(1) + qJ(3);
	t46 = cos(pkin(8));
	t45 = sin(pkin(8));
	t1 = [t48 * t45, t48 * t46, t49, qJ(2) * t48 + t47 * t49 + 0; -t49 * t45, -t49 * t46, t48, -qJ(2) * t49 + t47 * t48 + 0; t46, -t45, 0, pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:40
	% EndTime: 2020-11-04 20:16:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t54 = pkin(8) + qJ(4);
	t53 = pkin(1) + pkin(6) + qJ(3);
	t52 = cos(t54);
	t51 = sin(t54);
	t50 = sin(pkin(8)) * pkin(3) + qJ(2);
	t1 = [t55 * t51, t55 * t52, t56, t50 * t55 + t53 * t56 + 0; -t56 * t51, -t56 * t52, t55, -t50 * t56 + t53 * t55 + 0; t52, -t51, 0, cos(pkin(8)) * pkin(3) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:40
	% EndTime: 2020-11-04 20:16:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (40->20), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t62 = sin(qJ(5));
	t63 = sin(qJ(1));
	t70 = t63 * t62;
	t64 = cos(qJ(5));
	t69 = t63 * t64;
	t65 = cos(qJ(1));
	t68 = t65 * t62;
	t67 = t65 * t64;
	t61 = pkin(8) + qJ(4);
	t58 = sin(t61);
	t59 = cos(t61);
	t66 = pkin(4) * t58 - pkin(7) * t59 + sin(pkin(8)) * pkin(3) + qJ(2);
	t60 = pkin(1) + pkin(6) + qJ(3);
	t1 = [t58 * t69 + t68, -t58 * t70 + t67, -t63 * t59, t60 * t65 + t66 * t63 + 0; -t58 * t67 + t70, t58 * t68 + t69, t65 * t59, t60 * t63 - t66 * t65 + 0; t59 * t64, -t59 * t62, t58, t59 * pkin(4) + t58 * pkin(7) + cos(pkin(8)) * pkin(3) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end