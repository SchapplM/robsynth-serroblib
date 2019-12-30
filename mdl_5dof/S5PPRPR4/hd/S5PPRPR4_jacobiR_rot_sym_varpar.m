% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRPR4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PPRPR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:17:49
	% EndTime: 2019-12-29 15:17:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:17:49
	% EndTime: 2019-12-29 15:17:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:17:49
	% EndTime: 2019-12-29 15:17:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:17:49
	% EndTime: 2019-12-29 15:17:49
	% DurationCPUTime: 0.03s
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
	% StartTime: 2019-12-29 15:17:54
	% EndTime: 2019-12-29 15:17:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (8->5), mult. (16->8), div. (0->0), fcn. (28->6), ass. (0->9)
	t21 = cos(qJ(3));
	t20 = sin(qJ(3));
	t19 = cos(pkin(7));
	t18 = cos(pkin(8));
	t17 = sin(pkin(7));
	t16 = sin(pkin(8));
	t15 = t17 * t21 - t19 * t20;
	t14 = -t17 * t20 - t19 * t21;
	t1 = [0, 0, t15 * t18, 0, 0; 0, 0, t14 * t18, 0, 0; 0, 0, 0, 0, 0; 0, 0, -t15 * t16, 0, 0; 0, 0, -t14 * t16, 0, 0; 0, 0, 0, 0, 0; 0, 0, -t14, 0, 0; 0, 0, t15, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:17:54
	% EndTime: 2019-12-29 15:17:54
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (27->9), mult. (28->8), div. (0->0), fcn. (50->6), ass. (0->14)
	t37 = cos(qJ(3));
	t36 = sin(qJ(3));
	t30 = sin(pkin(7));
	t31 = cos(pkin(7));
	t23 = -t30 * t36 - t31 * t37;
	t29 = pkin(8) + qJ(5);
	t27 = sin(t29);
	t35 = t23 * t27;
	t28 = cos(t29);
	t34 = t23 * t28;
	t24 = t30 * t37 - t31 * t36;
	t33 = t24 * t27;
	t32 = t24 * t28;
	t1 = [0, 0, t32, 0, t35; 0, 0, t34, 0, -t33; 0, 0, 0, 0, -t28; 0, 0, -t33, 0, t34; 0, 0, -t35, 0, -t32; 0, 0, 0, 0, t27; 0, 0, -t23, 0, 0; 0, 0, t24, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end