% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRRP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:49
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:38
	% EndTime: 2020-11-04 19:49:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:38
	% EndTime: 2020-11-04 19:49:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = cos(qJ(1));
	t30 = sin(qJ(1));
	t1 = [t31, -t30, 0, 0; t30, t31, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:38
	% EndTime: 2020-11-04 19:49:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t35 = cos(qJ(1));
	t34 = cos(qJ(2));
	t33 = sin(qJ(1));
	t32 = sin(qJ(2));
	t1 = [t35 * t34, -t35 * t32, t33, t35 * pkin(1) + t33 * pkin(5) + 0; t33 * t34, -t33 * t32, -t35, t33 * pkin(1) - t35 * pkin(5) + 0; t32, t34, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:38
	% EndTime: 2020-11-04 19:49:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t42 = pkin(6) + pkin(5);
	t41 = cos(qJ(1));
	t40 = sin(qJ(1));
	t39 = qJ(2) + qJ(3);
	t38 = cos(t39);
	t37 = sin(t39);
	t36 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t41 * t38, -t41 * t37, t40, t41 * t36 + t42 * t40 + 0; t40 * t38, -t40 * t37, -t41, t40 * t36 - t41 * t42 + 0; t37, t38, 0, sin(qJ(2)) * pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:38
	% EndTime: 2020-11-04 19:49:38
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (30->15), mult. (21->14), div. (0->0), fcn. (29->6), ass. (0->8)
	t46 = qJ(2) + qJ(3);
	t44 = sin(t46);
	t45 = cos(t46);
	t50 = pkin(3) * t45 + qJ(4) * t44 + cos(qJ(2)) * pkin(2) + pkin(1);
	t49 = pkin(6) + pkin(5);
	t48 = cos(qJ(1));
	t47 = sin(qJ(1));
	t1 = [t48 * t45, t47, t48 * t44, t49 * t47 + t50 * t48 + 0; t47 * t45, -t48, t47 * t44, t50 * t47 - t48 * t49 + 0; t44, 0, -t45, t44 * pkin(3) - t45 * qJ(4) + sin(qJ(2)) * pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end