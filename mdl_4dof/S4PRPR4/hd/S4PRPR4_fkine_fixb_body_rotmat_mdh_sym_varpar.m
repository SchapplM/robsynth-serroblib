% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRPR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:35
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:35:02
	% EndTime: 2020-11-04 19:35:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:35:02
	% EndTime: 2020-11-04 19:35:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t25 = cos(pkin(6));
	t24 = sin(pkin(6));
	t1 = [t25, -t24, 0, 0; t24, t25, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:35:02
	% EndTime: 2020-11-04 19:35:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t28 = pkin(6) + qJ(2);
	t27 = cos(t28);
	t26 = sin(t28);
	t1 = [t27, -t26, 0, cos(pkin(6)) * pkin(1) + 0; t26, t27, 0, sin(pkin(6)) * pkin(1) + 0; 0, 0, 1, pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:35:02
	% EndTime: 2020-11-04 19:35:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (6->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t31 = pkin(6) + qJ(2);
	t30 = cos(t31);
	t29 = sin(t31);
	t1 = [0, -t30, t29, t30 * pkin(2) + t29 * qJ(3) + cos(pkin(6)) * pkin(1) + 0; 0, -t29, -t30, t29 * pkin(2) - t30 * qJ(3) + sin(pkin(6)) * pkin(1) + 0; 1, 0, 0, pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:35:02
	% EndTime: 2020-11-04 19:35:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (24->14), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->7)
	t37 = pkin(2) + pkin(5);
	t36 = cos(qJ(4));
	t35 = sin(qJ(4));
	t34 = pkin(6) + qJ(2);
	t33 = cos(t34);
	t32 = sin(t34);
	t1 = [t32 * t35, t32 * t36, t33, t37 * t33 + t32 * qJ(3) + cos(pkin(6)) * pkin(1) + 0; -t33 * t35, -t33 * t36, t32, t37 * t32 - t33 * qJ(3) + sin(pkin(6)) * pkin(1) + 0; t36, -t35, 0, pkin(3) + pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end