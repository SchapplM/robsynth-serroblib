% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:40
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:40:45
	% EndTime: 2020-11-04 20:40:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:40:45
	% EndTime: 2020-11-04 20:40:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t37 = cos(qJ(1));
	t36 = sin(qJ(1));
	t1 = [t37, -t36, 0, 0; t36, t37, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:40:45
	% EndTime: 2020-11-04 20:40:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t41 = cos(qJ(1));
	t40 = cos(qJ(2));
	t39 = sin(qJ(1));
	t38 = sin(qJ(2));
	t1 = [t41 * t40, -t41 * t38, t39, t41 * pkin(1) + t39 * pkin(6) + 0; t39 * t40, -t39 * t38, -t41, t39 * pkin(1) - t41 * pkin(6) + 0; t38, t40, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:40:45
	% EndTime: 2020-11-04 20:40:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t48 = pkin(7) + pkin(6);
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t45 = qJ(2) + qJ(3);
	t44 = cos(t45);
	t43 = sin(t45);
	t42 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t47 * t44, -t47 * t43, t46, t47 * t42 + t48 * t46 + 0; t46 * t44, -t46 * t43, -t47, t46 * t42 - t47 * t48 + 0; t43, t44, 0, sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:40:45
	% EndTime: 2020-11-04 20:40:45
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (30->15), mult. (21->14), div. (0->0), fcn. (29->6), ass. (0->8)
	t52 = qJ(2) + qJ(3);
	t50 = sin(t52);
	t51 = cos(t52);
	t56 = pkin(3) * t51 + qJ(4) * t50 + cos(qJ(2)) * pkin(2) + pkin(1);
	t55 = pkin(7) + pkin(6);
	t54 = cos(qJ(1));
	t53 = sin(qJ(1));
	t1 = [t54 * t51, t53, t54 * t50, t55 * t53 + t56 * t54 + 0; t53 * t51, -t54, t53 * t50, t56 * t53 - t54 * t55 + 0; t50, 0, -t51, t50 * pkin(3) - t51 * qJ(4) + sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:40:45
	% EndTime: 2020-11-04 20:40:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (37->18), mult. (25->17), div. (0->0), fcn. (33->8), ass. (0->12)
	t67 = pkin(3) + pkin(4);
	t66 = cos(qJ(1));
	t65 = cos(qJ(3));
	t64 = sin(qJ(1));
	t63 = sin(qJ(2));
	t62 = sin(qJ(3));
	t61 = qJ(2) + qJ(3);
	t60 = qJ(5) - pkin(6) - pkin(7);
	t59 = cos(t61);
	t58 = sin(t61);
	t57 = (qJ(4) * t62 + t67 * t65 + pkin(2)) * cos(qJ(2)) + pkin(1) + (qJ(4) * t65 - t62 * t67) * t63;
	t1 = [t66 * t59, t66 * t58, -t64, t57 * t66 - t60 * t64 + 0; t64 * t59, t64 * t58, t66, t57 * t64 + t60 * t66 + 0; t58, -t59, 0, t63 * pkin(2) - t59 * qJ(4) + t67 * t58 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end