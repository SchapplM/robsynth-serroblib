% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRPP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:46
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:46:07
	% EndTime: 2020-11-04 19:46:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:46:07
	% EndTime: 2020-11-04 19:46:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t25 = cos(qJ(1));
	t24 = sin(qJ(1));
	t1 = [t25, -t24, 0, 0; t24, t25, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:46:07
	% EndTime: 2020-11-04 19:46:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t29 = cos(qJ(1));
	t28 = cos(qJ(2));
	t27 = sin(qJ(1));
	t26 = sin(qJ(2));
	t1 = [t29 * t28, -t29 * t26, t27, t29 * pkin(1) + t27 * pkin(5) + 0; t27 * t28, -t27 * t26, -t29, t27 * pkin(1) - t29 * pkin(5) + 0; t26, t28, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:46:07
	% EndTime: 2020-11-04 19:46:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->14), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t34 = cos(qJ(1));
	t33 = cos(qJ(2));
	t32 = sin(qJ(1));
	t31 = sin(qJ(2));
	t30 = t33 * pkin(2) + t31 * qJ(3) + pkin(1);
	t1 = [t32, -t34 * t33, t34 * t31, t32 * pkin(5) + t30 * t34 + 0; -t34, -t32 * t33, t32 * t31, -t34 * pkin(5) + t30 * t32 + 0; 0, -t31, -t33, t31 * pkin(2) - t33 * qJ(3) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:46:07
	% EndTime: 2020-11-04 19:46:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->13), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->8)
	t41 = pkin(3) + pkin(5);
	t40 = cos(qJ(1));
	t39 = cos(qJ(2));
	t38 = sin(qJ(1));
	t37 = sin(qJ(2));
	t36 = pkin(2) + qJ(4);
	t35 = t37 * qJ(3) + t36 * t39 + pkin(1);
	t1 = [t38, t40 * t37, t40 * t39, t35 * t40 + t41 * t38 + 0; -t40, t38 * t37, t38 * t39, t35 * t38 - t41 * t40 + 0; 0, -t39, t37, -t39 * qJ(3) + t36 * t37 + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end