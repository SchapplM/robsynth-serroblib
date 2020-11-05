% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRPR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:48
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:47:59
	% EndTime: 2020-11-04 19:47:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:47:59
	% EndTime: 2020-11-04 19:47:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t33 = cos(qJ(1));
	t32 = sin(qJ(1));
	t1 = [t33, -t32, 0, 0; t32, t33, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:47:59
	% EndTime: 2020-11-04 19:47:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t37 = cos(qJ(1));
	t36 = cos(qJ(2));
	t35 = sin(qJ(1));
	t34 = sin(qJ(2));
	t1 = [t37 * t36, -t37 * t34, t35, t37 * pkin(1) + t35 * pkin(5) + 0; t35 * t36, -t35 * t34, -t37, t35 * pkin(1) - t37 * pkin(5) + 0; t34, t36, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:47:59
	% EndTime: 2020-11-04 19:47:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t42 = cos(qJ(1));
	t41 = cos(qJ(2));
	t40 = sin(qJ(1));
	t39 = sin(qJ(2));
	t38 = t41 * pkin(2) + t39 * qJ(3) + pkin(1);
	t1 = [t42 * t41, t40, t42 * t39, t40 * pkin(5) + t38 * t42 + 0; t40 * t41, -t42, t40 * t39, -t42 * pkin(5) + t38 * t40 + 0; t39, 0, -t41, t39 * pkin(2) - t41 * qJ(3) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:47:59
	% EndTime: 2020-11-04 19:47:59
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (23->16), mult. (26->16), div. (0->0), fcn. (40->6), ass. (0->12)
	t55 = cos(qJ(4));
	t54 = sin(qJ(4));
	t53 = pkin(2) + pkin(3);
	t52 = pkin(5) - pkin(6);
	t51 = cos(qJ(1));
	t50 = cos(qJ(2));
	t49 = sin(qJ(1));
	t48 = sin(qJ(2));
	t45 = t48 * qJ(3) + t53 * t50 + pkin(1);
	t44 = t48 * t55 - t50 * t54;
	t43 = -t48 * t54 - t50 * t55;
	t1 = [-t51 * t43, t51 * t44, -t49, t45 * t51 + t52 * t49 + 0; -t49 * t43, t49 * t44, t51, t45 * t49 - t52 * t51 + 0; t44, t43, 0, -t50 * qJ(3) + t53 * t48 + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end