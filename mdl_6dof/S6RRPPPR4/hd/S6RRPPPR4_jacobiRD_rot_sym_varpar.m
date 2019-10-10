% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:23
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPPR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t34 = sin(qJ(1));
	t41 = qJD(1) * t34;
	t36 = cos(qJ(1));
	t40 = qJD(1) * t36;
	t33 = sin(qJ(2));
	t39 = qJD(2) * t33;
	t35 = cos(qJ(2));
	t38 = qJD(2) * t35;
	t37 = qJD(2) * t36;
	t32 = t34 * t39 - t35 * t40;
	t31 = t33 * t40 + t34 * t38;
	t30 = t33 * t37 + t35 * t41;
	t29 = t33 * t41 - t35 * t37;
	t1 = [t32, t29, 0, 0, 0, 0; -t30, -t31, 0, 0, 0, 0; 0, -t39, 0, 0, 0, 0; t31, t30, 0, 0, 0, 0; t29, t32, 0, 0, 0, 0; 0, -t38, 0, 0, 0, 0; -t41, 0, 0, 0, 0, 0; t40, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:25
	% EndTime: 2019-10-10 09:23:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t161 = sin(qJ(1));
	t168 = qJD(1) * t161;
	t163 = cos(qJ(1));
	t167 = qJD(1) * t163;
	t160 = sin(qJ(2));
	t166 = qJD(2) * t160;
	t162 = cos(qJ(2));
	t165 = qJD(2) * t162;
	t164 = qJD(2) * t163;
	t159 = -t161 * t166 + t162 * t167;
	t158 = t160 * t167 + t161 * t165;
	t157 = t160 * t164 + t162 * t168;
	t156 = -t160 * t168 + t162 * t164;
	t1 = [-t168, 0, 0, 0, 0, 0; t167, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t159, t156, 0, 0, 0, 0; t157, t158, 0, 0, 0, 0; 0, t166, 0, 0, 0, 0; -t158, -t157, 0, 0, 0, 0; t156, t159, 0, 0, 0, 0; 0, t165, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:25
	% EndTime: 2019-10-10 09:23:25
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (17->14), mult. (77->37), div. (0->0), fcn. (77->6), ass. (0->18)
	t187 = sin(qJ(2));
	t188 = sin(qJ(1));
	t201 = t187 * t188;
	t190 = cos(qJ(1));
	t200 = t187 * t190;
	t199 = qJD(1) * t188;
	t198 = qJD(1) * t190;
	t197 = qJD(2) * t187;
	t189 = cos(qJ(2));
	t196 = qJD(2) * t189;
	t195 = qJD(2) * t190;
	t194 = t188 * t196;
	t193 = t189 * t195;
	t192 = -t188 * t197 + t189 * t198;
	t191 = -t187 * t195 - t189 * t199;
	t186 = cos(pkin(9));
	t185 = sin(pkin(9));
	t1 = [-t185 * t194 + (-t185 * t200 - t186 * t188) * qJD(1), t191 * t185, 0, 0, 0, 0; t185 * t193 + (-t185 * t201 + t186 * t190) * qJD(1), t192 * t185, 0, 0, 0, 0; 0, t185 * t196, 0, 0, 0, 0; -t186 * t194 + (t185 * t188 - t186 * t200) * qJD(1), t191 * t186, 0, 0, 0, 0; t186 * t193 + (-t185 * t190 - t186 * t201) * qJD(1), t192 * t186, 0, 0, 0, 0; 0, t186 * t196, 0, 0, 0, 0; -t192, t187 * t199 - t193, 0, 0, 0, 0; t191, -t187 * t198 - t194, 0, 0, 0, 0; 0, -t197, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:25
	% EndTime: 2019-10-10 09:23:25
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (18->18), mult. (77->37), div. (0->0), fcn. (77->6), ass. (0->18)
	t209 = sin(qJ(2));
	t210 = sin(qJ(1));
	t223 = t209 * t210;
	t212 = cos(qJ(1));
	t222 = t209 * t212;
	t221 = qJD(1) * t210;
	t220 = qJD(1) * t212;
	t219 = qJD(2) * t209;
	t211 = cos(qJ(2));
	t218 = qJD(2) * t211;
	t217 = qJD(2) * t212;
	t216 = t210 * t218;
	t215 = t211 * t217;
	t214 = -t210 * t219 + t211 * t220;
	t213 = t209 * t217 + t211 * t221;
	t208 = cos(pkin(9));
	t207 = sin(pkin(9));
	t1 = [-t207 * t216 + (-t207 * t222 - t208 * t210) * qJD(1), -t213 * t207, 0, 0, 0, 0; t207 * t215 + (-t207 * t223 + t208 * t212) * qJD(1), t214 * t207, 0, 0, 0, 0; 0, t207 * t218, 0, 0, 0, 0; -t214, t209 * t221 - t215, 0, 0, 0, 0; -t213, -t209 * t220 - t216, 0, 0, 0, 0; 0, -t219, 0, 0, 0, 0; t208 * t216 + (-t207 * t210 + t208 * t222) * qJD(1), t213 * t208, 0, 0, 0, 0; -t208 * t215 + (t207 * t212 + t208 * t223) * qJD(1), -t214 * t208, 0, 0, 0, 0; 0, -t208 * t218, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:26
	% EndTime: 2019-10-10 09:23:26
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (108->39), mult. (385->79), div. (0->0), fcn. (401->8), ass. (0->43)
	t305 = sin(qJ(2));
	t308 = cos(qJ(2));
	t306 = sin(qJ(1));
	t327 = qJD(2) * t306;
	t323 = t308 * t327;
	t309 = cos(qJ(1));
	t329 = qJD(1) * t309;
	t336 = t305 * t329 + t323;
	t302 = sin(pkin(9));
	t303 = cos(pkin(9));
	t330 = qJD(1) * t306;
	t291 = t302 * t330 - t336 * t303;
	t332 = t309 * t302;
	t333 = t306 * t303;
	t296 = t305 * t332 + t333;
	t292 = t296 * qJD(1) + t302 * t323;
	t297 = t305 * t333 + t332;
	t331 = t309 * t303;
	t334 = t306 * t302;
	t298 = -t305 * t334 + t331;
	t304 = sin(qJ(6));
	t307 = cos(qJ(6));
	t335 = t291 * t307 - t292 * t304 + (t297 * t304 + t298 * t307) * qJD(6);
	t328 = qJD(2) * t305;
	t326 = qJD(2) * t309;
	t324 = t308 * t329;
	t322 = t308 * t326;
	t319 = t302 * t307 - t303 * t304;
	t318 = t302 * t304 + t303 * t307;
	t317 = t319 * t308;
	t316 = qJD(2) * t318;
	t315 = qJD(6) * t319;
	t314 = qJD(6) * t318;
	t313 = t319 * t328;
	t312 = t318 * t328;
	t311 = -t291 * t304 - t292 * t307 + (t297 * t307 - t298 * t304) * qJD(6);
	t310 = -t308 * t314 - t313;
	t295 = -t305 * t331 + t334;
	t294 = t298 * qJD(1) + t302 * t322;
	t293 = t297 * qJD(1) - t303 * t322;
	t290 = t293 * t304 + t294 * t307 + (t295 * t307 - t296 * t304) * qJD(6);
	t289 = t293 * t307 - t294 * t304 + (-t295 * t304 - t296 * t307) * qJD(6);
	t1 = [t311, -t309 * t313 + (-t309 * t314 - t319 * t330) * t308, 0, 0, 0, t289; t290, t310 * t306 + t317 * t329, 0, 0, 0, t335; 0, qJD(2) * t317 - t305 * t314, 0, 0, 0, qJD(6) * t317 - t305 * t316; -t335, t309 * t312 + (-t309 * t315 + t318 * t330) * t308, 0, 0, 0, -t290; t289, -t318 * t324 + (-t308 * t315 + t312) * t306, 0, 0, 0, t311; 0, -t305 * t315 - t308 * t316, 0, 0, 0, t310; -t305 * t327 + t324, -t305 * t330 + t322, 0, 0, 0, 0; t305 * t326 + t308 * t330, t336, 0, 0, 0, 0; 0, t328, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end