% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:17
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:28
	% EndTime: 2020-11-04 20:17:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:28
	% EndTime: 2020-11-04 20:17:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t36 = cos(qJ(1));
	t35 = sin(qJ(1));
	t1 = [t36, -t35, 0, 0; t35, t36, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:28
	% EndTime: 2020-11-04 20:17:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t40 = cos(qJ(1));
	t39 = sin(qJ(1));
	t38 = cos(pkin(7));
	t37 = sin(pkin(7));
	t1 = [t40 * t38, -t40 * t37, t39, t40 * pkin(1) + t39 * qJ(2) + 0; t39 * t38, -t39 * t37, -t40, t39 * pkin(1) - t40 * qJ(2) + 0; t37, t38, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:28
	% EndTime: 2020-11-04 20:17:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t45 = pkin(6) + qJ(2);
	t44 = pkin(7) + qJ(3);
	t43 = cos(t44);
	t42 = sin(t44);
	t41 = cos(pkin(7)) * pkin(2) + pkin(1);
	t1 = [t47 * t43, -t47 * t42, t46, t47 * t41 + t45 * t46 + 0; t46 * t43, -t46 * t42, -t47, t46 * t41 - t47 * t45 + 0; t42, t43, 0, sin(pkin(7)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:28
	% EndTime: 2020-11-04 20:17:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->18), mult. (21->14), div. (0->0), fcn. (29->6), ass. (0->8)
	t51 = pkin(7) + qJ(3);
	t49 = sin(t51);
	t50 = cos(t51);
	t55 = pkin(3) * t50 + qJ(4) * t49 + cos(pkin(7)) * pkin(2) + pkin(1);
	t54 = cos(qJ(1));
	t53 = sin(qJ(1));
	t52 = pkin(6) + qJ(2);
	t1 = [t53, -t54 * t50, t54 * t49, t52 * t53 + t55 * t54 + 0; -t54, -t53 * t50, t53 * t49, -t54 * t52 + t55 * t53 + 0; 0, -t49, -t50, t49 * pkin(3) - t50 * qJ(4) + sin(pkin(7)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:28
	% EndTime: 2020-11-04 20:17:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (37->18), mult. (25->18), div. (0->0), fcn. (33->8), ass. (0->11)
	t65 = cos(qJ(1));
	t64 = sin(qJ(1));
	t63 = pkin(3) + qJ(5);
	t62 = cos(pkin(7));
	t61 = sin(pkin(7));
	t60 = pkin(7) + qJ(3);
	t59 = pkin(4) + pkin(6) + qJ(2);
	t58 = cos(t60);
	t57 = sin(t60);
	t56 = (qJ(4) * t61 + t63 * t62) * cos(qJ(3)) + (qJ(4) * t62 - t61 * t63) * sin(qJ(3)) + t62 * pkin(2) + pkin(1);
	t1 = [t64, t65 * t57, t65 * t58, t56 * t65 + t59 * t64 + 0; -t65, t64 * t57, t64 * t58, t56 * t64 - t59 * t65 + 0; 0, -t58, t57, t61 * pkin(2) - t58 * qJ(4) + t63 * t57 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end