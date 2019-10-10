% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPPRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(9);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0, 0; -t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t37, 0, 0, 0, 0, 0; -t36, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (5->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t16 = qJ(1) + pkin(9);
	t17 = qJD(1) * sin(t16);
	t14 = qJD(1) * cos(t16);
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t14, 0, 0, 0, 0, 0; t17, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t17, 0, 0, 0, 0, 0; t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t18 = qJ(1) + pkin(9);
	t20 = qJD(1) * sin(t18);
	t19 = qJD(1) * cos(t18);
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t20, 0, 0, 0, 0, 0; t19, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t19, 0, 0, 0, 0, 0; -t20, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:53
	% EndTime: 2019-10-09 23:25:53
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (28->9), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t44 = sin(qJ(5));
	t49 = qJD(1) * t44;
	t45 = cos(qJ(5));
	t48 = qJD(1) * t45;
	t47 = qJD(5) * t44;
	t46 = qJD(5) * t45;
	t43 = qJ(1) + pkin(9);
	t42 = cos(t43);
	t41 = sin(t43);
	t40 = -t41 * t47 + t42 * t48;
	t39 = -t41 * t46 - t42 * t49;
	t38 = -t41 * t48 - t42 * t47;
	t37 = t41 * t49 - t42 * t46;
	t1 = [t39, 0, 0, 0, t38, 0; -t37, 0, 0, 0, t40, 0; 0, 0, 0, 0, -t46, 0; -t40, 0, 0, 0, t37, 0; t38, 0, 0, 0, t39, 0; 0, 0, 0, 0, t47, 0; qJD(1) * t41, 0, 0, 0, 0, 0; -qJD(1) * t42, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:54
	% EndTime: 2019-10-09 23:25:54
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (109->24), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->30)
	t244 = cos(qJ(6));
	t243 = sin(qJ(5));
	t261 = qJD(1) * t243;
	t250 = qJD(6) + t261;
	t263 = t244 * t250;
	t257 = qJD(6) * t243;
	t251 = qJD(1) + t257;
	t242 = sin(qJ(6));
	t245 = cos(qJ(5));
	t258 = qJD(5) * t245;
	t253 = t242 * t258;
	t262 = t251 * t244 + t253;
	t260 = qJD(1) * t245;
	t259 = qJD(5) * t243;
	t256 = qJD(6) * t245;
	t255 = t242 * t260;
	t254 = t244 * t260;
	t252 = t244 * t258;
	t249 = t250 * t242;
	t248 = t242 * t256 + t244 * t259;
	t247 = t242 * t259 - t244 * t256;
	t246 = t251 * t242 - t252;
	t241 = qJ(1) + pkin(9);
	t240 = cos(t241);
	t239 = sin(t241);
	t238 = t246 * t239 - t240 * t263;
	t237 = t262 * t239 + t240 * t249;
	t236 = t239 * t263 + t246 * t240;
	t235 = t239 * t249 - t262 * t240;
	t1 = [t238, 0, 0, 0, -t239 * t254 - t248 * t240, t235; -t236, 0, 0, 0, -t248 * t239 + t240 * t254, -t237; 0, 0, 0, 0, t242 * t257 - t252, t247; t237, 0, 0, 0, t239 * t255 + t247 * t240, t236; t235, 0, 0, 0, t247 * t239 - t240 * t255, t238; 0, 0, 0, 0, t244 * t257 + t253, t248; -t239 * t259 + t240 * t260, 0, 0, 0, -t239 * t261 + t240 * t258, 0; t239 * t260 + t240 * t259, 0, 0, 0, t239 * t258 + t240 * t261, 0; 0, 0, 0, 0, -t259, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end