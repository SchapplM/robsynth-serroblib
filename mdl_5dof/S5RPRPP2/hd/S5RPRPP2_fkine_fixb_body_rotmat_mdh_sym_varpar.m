% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPP2 (for one body)
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

function Tc_mdh = S5RPRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:12
	% EndTime: 2020-11-04 20:17:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:12
	% EndTime: 2020-11-04 20:17:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t40 = cos(qJ(1));
	t39 = sin(qJ(1));
	t1 = [t40, -t39, 0, 0; t39, t40, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:12
	% EndTime: 2020-11-04 20:17:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t43 = qJ(1) + pkin(7);
	t42 = cos(t43);
	t41 = sin(t43);
	t1 = [t42, -t41, 0, cos(qJ(1)) * pkin(1) + 0; t41, t42, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:12
	% EndTime: 2020-11-04 20:17:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t48 = cos(qJ(3));
	t47 = sin(qJ(3));
	t46 = qJ(1) + pkin(7);
	t45 = cos(t46);
	t44 = sin(t46);
	t1 = [t45 * t48, -t45 * t47, t44, t45 * pkin(2) + t44 * pkin(6) + cos(qJ(1)) * pkin(1) + 0; t44 * t48, -t44 * t47, -t45, t44 * pkin(2) - t45 * pkin(6) + sin(qJ(1)) * pkin(1) + 0; t47, t48, 0, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:12
	% EndTime: 2020-11-04 20:17:12
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (30->15), mult. (20->14), div. (0->0), fcn. (28->6), ass. (0->7)
	t52 = sin(qJ(3));
	t53 = cos(qJ(3));
	t54 = pkin(3) * t53 + qJ(4) * t52 + pkin(2);
	t51 = qJ(1) + pkin(7);
	t50 = cos(t51);
	t49 = sin(t51);
	t1 = [t50 * t53, t49, t50 * t52, cos(qJ(1)) * pkin(1) + t49 * pkin(6) + 0 + t54 * t50; t49 * t53, -t50, t49 * t52, sin(qJ(1)) * pkin(1) - t50 * pkin(6) + 0 + t54 * t49; t52, 0, -t53, t52 * pkin(3) - qJ(4) * t53 + pkin(5) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:12
	% EndTime: 2020-11-04 20:17:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (39->17), mult. (27->14), div. (0->0), fcn. (35->6), ass. (0->9)
	t62 = pkin(3) + pkin(4);
	t61 = pkin(6) - qJ(5);
	t58 = sin(qJ(3));
	t59 = cos(qJ(3));
	t60 = qJ(4) * t58 + t62 * t59 + pkin(2);
	t57 = qJ(1) + pkin(7);
	t56 = cos(t57);
	t55 = sin(t57);
	t1 = [t56 * t59, t56 * t58, -t55, cos(qJ(1)) * pkin(1) + 0 + t61 * t55 + t60 * t56; t55 * t59, t55 * t58, t56, sin(qJ(1)) * pkin(1) + 0 - t61 * t56 + t60 * t55; t58, -t59, 0, -t59 * qJ(4) + t62 * t58 + pkin(5) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end