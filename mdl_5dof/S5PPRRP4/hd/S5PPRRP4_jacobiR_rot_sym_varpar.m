% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PPRRP4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_jacobiR_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:34:52
	% EndTime: 2019-12-31 17:34:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:34:52
	% EndTime: 2019-12-31 17:34:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:34:52
	% EndTime: 2019-12-31 17:34:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:34:52
	% EndTime: 2019-12-31 17:34:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->3), mult. (8->4), div. (0->0), fcn. (16->4), ass. (0->7)
	t20 = cos(qJ(3));
	t19 = sin(qJ(3));
	t18 = cos(pkin(7));
	t17 = sin(pkin(7));
	t16 = -t17 * t20 + t18 * t19;
	t15 = -t17 * t19 - t18 * t20;
	t1 = [0, 0, -t16, 0, 0; 0, 0, t15, 0, 0; 0, 0, 0, 0, 0; 0, 0, t15, 0, 0; 0, 0, t16, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:34:52
	% EndTime: 2019-12-31 17:34:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (28->8), div. (0->0), fcn. (50->6), ass. (0->13)
	t33 = cos(qJ(3));
	t32 = sin(qJ(3));
	t26 = sin(pkin(7));
	t27 = cos(pkin(7));
	t20 = -t26 * t32 - t27 * t33;
	t24 = sin(qJ(4));
	t31 = t20 * t24;
	t25 = cos(qJ(4));
	t30 = t20 * t25;
	t21 = t26 * t33 - t27 * t32;
	t29 = t21 * t24;
	t28 = t21 * t25;
	t1 = [0, 0, t28, t31, 0; 0, 0, t30, -t29, 0; 0, 0, 0, -t25, 0; 0, 0, -t29, t30, 0; 0, 0, -t31, -t28, 0; 0, 0, 0, t24, 0; 0, 0, -t20, 0, 0; 0, 0, t21, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:34:52
	% EndTime: 2019-12-31 17:34:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (28->8), div. (0->0), fcn. (50->6), ass. (0->13)
	t37 = cos(qJ(3));
	t36 = sin(qJ(3));
	t30 = sin(pkin(7));
	t31 = cos(pkin(7));
	t24 = -t30 * t36 - t31 * t37;
	t28 = sin(qJ(4));
	t35 = t24 * t28;
	t29 = cos(qJ(4));
	t34 = t24 * t29;
	t25 = t30 * t37 - t31 * t36;
	t33 = t25 * t28;
	t32 = t25 * t29;
	t1 = [0, 0, t32, t35, 0; 0, 0, t34, -t33, 0; 0, 0, 0, -t29, 0; 0, 0, -t33, t34, 0; 0, 0, -t35, -t32, 0; 0, 0, 0, t28, 0; 0, 0, -t24, 0, 0; 0, 0, t25, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end