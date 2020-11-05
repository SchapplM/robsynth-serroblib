% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:33
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:33:35
	% EndTime: 2020-11-04 20:33:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:33:35
	% EndTime: 2020-11-04 20:33:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t40 = cos(qJ(1));
	t39 = sin(qJ(1));
	t1 = [t40, -t39, 0, 0; t39, t40, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:33:35
	% EndTime: 2020-11-04 20:33:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t44 = cos(qJ(1));
	t43 = cos(qJ(2));
	t42 = sin(qJ(1));
	t41 = sin(qJ(2));
	t1 = [t44 * t43, -t44 * t41, t42, t44 * pkin(1) + t42 * pkin(6) + 0; t42 * t43, -t42 * t41, -t44, t42 * pkin(1) - t44 * pkin(6) + 0; t41, t43, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:33:35
	% EndTime: 2020-11-04 20:33:35
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t49 = -qJ(3) - pkin(6);
	t48 = qJ(2) + pkin(8);
	t47 = cos(t48);
	t46 = sin(t48);
	t45 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t51 * t47, -t51 * t46, t50, t51 * t45 - t49 * t50 + 0; t50 * t47, -t50 * t46, -t51, t50 * t45 + t51 * t49 + 0; t46, t47, 0, sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:33:35
	% EndTime: 2020-11-04 20:33:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->15), mult. (14->12), div. (0->0), fcn. (22->8), ass. (0->9)
	t57 = qJ(2) + pkin(8);
	t59 = cos(qJ(1));
	t58 = sin(qJ(1));
	t56 = -pkin(7) - qJ(3) - pkin(6);
	t55 = qJ(4) + t57;
	t54 = cos(t55);
	t53 = sin(t55);
	t52 = pkin(3) * cos(t57) + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t59 * t54, -t59 * t53, t58, t59 * t52 - t58 * t56 + 0; t58 * t54, -t58 * t53, -t59, t58 * t52 + t59 * t56 + 0; t53, t54, 0, pkin(3) * sin(t57) + sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:33:35
	% EndTime: 2020-11-04 20:33:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (50->18), mult. (24->16), div. (0->0), fcn. (32->8), ass. (0->9)
	t65 = qJ(2) + pkin(8);
	t63 = qJ(4) + t65;
	t61 = sin(t63);
	t62 = cos(t63);
	t68 = pkin(4) * t62 + qJ(5) * t61 + pkin(3) * cos(t65) + cos(qJ(2)) * pkin(2) + pkin(1);
	t67 = cos(qJ(1));
	t66 = sin(qJ(1));
	t64 = -pkin(7) - qJ(3) - pkin(6);
	t1 = [t67 * t62, t66, t67 * t61, -t66 * t64 + t68 * t67 + 0; t66 * t62, -t67, t66 * t61, t67 * t64 + t68 * t66 + 0; t61, 0, -t62, t61 * pkin(4) - t62 * qJ(5) + pkin(3) * sin(t65) + sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end