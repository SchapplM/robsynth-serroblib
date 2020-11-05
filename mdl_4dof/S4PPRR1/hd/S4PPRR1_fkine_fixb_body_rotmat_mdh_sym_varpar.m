% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PPRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:32
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:32:01
	% EndTime: 2020-11-04 19:32:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:32:01
	% EndTime: 2020-11-04 19:32:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t30 = cos(pkin(6));
	t29 = sin(pkin(6));
	t1 = [t30, -t29, 0, 0; t29, t30, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:32:01
	% EndTime: 2020-11-04 19:32:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t32 = cos(pkin(6));
	t31 = sin(pkin(6));
	t1 = [t32, 0, t31, t32 * pkin(1) + t31 * qJ(2) + 0; t31, 0, -t32, t31 * pkin(1) - t32 * qJ(2) + 0; 0, 1, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:32:01
	% EndTime: 2020-11-04 19:32:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t39 = pkin(1) + pkin(2);
	t38 = cos(qJ(3));
	t37 = sin(qJ(3));
	t36 = cos(pkin(6));
	t35 = sin(pkin(6));
	t34 = t35 * t38 - t36 * t37;
	t33 = -t35 * t37 - t36 * t38;
	t1 = [-t33, t34, 0, t35 * qJ(2) + t39 * t36 + 0; t34, t33, 0, -t36 * qJ(2) + t39 * t35 + 0; 0, 0, -1, -pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:32:01
	% EndTime: 2020-11-04 19:32:01
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (25->16), mult. (20->14), div. (0->0), fcn. (28->6), ass. (0->11)
	t49 = pkin(1) + pkin(2);
	t48 = cos(qJ(3));
	t47 = sin(qJ(3));
	t46 = cos(pkin(6));
	t45 = sin(pkin(6));
	t44 = qJ(3) + qJ(4);
	t43 = cos(t44);
	t42 = sin(t44);
	t41 = -t46 * t42 + t45 * t43;
	t40 = -t45 * t42 - t46 * t43;
	t1 = [-t40, t41, 0, t45 * qJ(2) + t49 * t46 + 0 + (t45 * t47 + t46 * t48) * pkin(3); t41, t40, 0, -t46 * qJ(2) + t49 * t45 + 0 + (t45 * t48 - t46 * t47) * pkin(3); 0, 0, -1, -pkin(5) - pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end