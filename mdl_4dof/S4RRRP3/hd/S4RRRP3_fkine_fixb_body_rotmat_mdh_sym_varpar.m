% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRRP3 (for one body)
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

function Tc_mdh = S4RRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:10
	% EndTime: 2020-11-04 19:49:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:10
	% EndTime: 2020-11-04 19:49:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = cos(qJ(1));
	t30 = sin(qJ(1));
	t1 = [t31, -t30, 0, 0; t30, t31, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:10
	% EndTime: 2020-11-04 19:49:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t34 = qJ(1) + qJ(2);
	t33 = cos(t34);
	t32 = sin(t34);
	t1 = [t33, -t32, 0, cos(qJ(1)) * pkin(1) + 0; t32, t33, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(5) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:10
	% EndTime: 2020-11-04 19:49:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t39 = cos(qJ(3));
	t38 = sin(qJ(3));
	t37 = qJ(1) + qJ(2);
	t36 = cos(t37);
	t35 = sin(t37);
	t1 = [t36 * t39, -t36 * t38, t35, t36 * pkin(2) + t35 * pkin(6) + cos(qJ(1)) * pkin(1) + 0; t35 * t39, -t35 * t38, -t36, t35 * pkin(2) - t36 * pkin(6) + sin(qJ(1)) * pkin(1) + 0; t38, t39, 0, pkin(5) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:10
	% EndTime: 2020-11-04 19:49:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->15), mult. (20->14), div. (0->0), fcn. (28->6), ass. (0->7)
	t43 = sin(qJ(3));
	t44 = cos(qJ(3));
	t45 = pkin(3) * t44 + qJ(4) * t43 + pkin(2);
	t42 = qJ(1) + qJ(2);
	t41 = cos(t42);
	t40 = sin(t42);
	t1 = [t41 * t44, t40, t41 * t43, cos(qJ(1)) * pkin(1) + t40 * pkin(6) + 0 + t45 * t41; t40 * t44, -t41, t40 * t43, sin(qJ(1)) * pkin(1) - t41 * pkin(6) + 0 + t45 * t40; t43, 0, -t44, t43 * pkin(3) - t44 * qJ(4) + pkin(4) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end