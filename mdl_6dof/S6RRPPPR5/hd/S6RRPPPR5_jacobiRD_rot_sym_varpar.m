% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:25
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPPR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
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
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
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
	% StartTime: 2019-10-10 09:25:10
	% EndTime: 2019-10-10 09:25:10
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (17->17), mult. (77->37), div. (0->0), fcn. (77->6), ass. (0->18)
	t180 = sin(qJ(1));
	t181 = cos(qJ(2));
	t193 = t180 * t181;
	t182 = cos(qJ(1));
	t192 = t181 * t182;
	t191 = qJD(1) * t180;
	t190 = qJD(1) * t182;
	t179 = sin(qJ(2));
	t189 = qJD(2) * t179;
	t188 = qJD(2) * t181;
	t187 = qJD(2) * t182;
	t186 = t180 * t189;
	t185 = t179 * t187;
	t184 = t179 * t190 + t180 * t188;
	t183 = t179 * t191 - t181 * t187;
	t178 = cos(pkin(9));
	t177 = sin(pkin(9));
	t1 = [t178 * t186 + (-t177 * t180 - t178 * t192) * qJD(1), t183 * t178, 0, 0, 0, 0; -t178 * t185 + (t177 * t182 - t178 * t193) * qJD(1), -t184 * t178, 0, 0, 0, 0; 0, -t178 * t189, 0, 0, 0, 0; -t177 * t186 + (t177 * t192 - t178 * t180) * qJD(1), -t183 * t177, 0, 0, 0, 0; t177 * t185 + (t177 * t193 + t178 * t182) * qJD(1), t184 * t177, 0, 0, 0, 0; 0, t177 * t189, 0, 0, 0, 0; -t184, -t181 * t191 - t185, 0, 0, 0, 0; -t183, t181 * t190 - t186, 0, 0, 0, 0; 0, t188, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:10
	% EndTime: 2019-10-10 09:25:10
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (17->17), mult. (77->37), div. (0->0), fcn. (77->6), ass. (0->18)
	t220 = sin(qJ(1));
	t221 = cos(qJ(2));
	t233 = t220 * t221;
	t222 = cos(qJ(1));
	t232 = t221 * t222;
	t231 = qJD(1) * t220;
	t230 = qJD(1) * t222;
	t219 = sin(qJ(2));
	t229 = qJD(2) * t219;
	t228 = qJD(2) * t221;
	t227 = qJD(2) * t222;
	t226 = t220 * t229;
	t225 = t219 * t227;
	t224 = t219 * t230 + t220 * t228;
	t223 = t219 * t231 - t221 * t227;
	t218 = cos(pkin(9));
	t217 = sin(pkin(9));
	t1 = [-t224, -t221 * t231 - t225, 0, 0, 0, 0; -t223, t221 * t230 - t226, 0, 0, 0, 0; 0, t228, 0, 0, 0, 0; -t218 * t226 + (t217 * t220 + t218 * t232) * qJD(1), -t223 * t218, 0, 0, 0, 0; t218 * t225 + (-t217 * t222 + t218 * t233) * qJD(1), t224 * t218, 0, 0, 0, 0; 0, t218 * t229, 0, 0, 0, 0; t217 * t226 + (-t217 * t232 + t218 * t220) * qJD(1), t223 * t217, 0, 0, 0, 0; -t217 * t225 + (-t217 * t233 - t218 * t222) * qJD(1), -t224 * t217, 0, 0, 0, 0; 0, -t217 * t229, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:10
	% EndTime: 2019-10-10 09:25:10
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (19->17), mult. (77->37), div. (0->0), fcn. (77->6), ass. (0->18)
	t213 = sin(qJ(1));
	t214 = cos(qJ(2));
	t226 = t213 * t214;
	t215 = cos(qJ(1));
	t225 = t214 * t215;
	t224 = qJD(1) * t213;
	t223 = qJD(1) * t215;
	t212 = sin(qJ(2));
	t222 = qJD(2) * t212;
	t221 = qJD(2) * t214;
	t220 = qJD(2) * t215;
	t219 = t213 * t222;
	t218 = t212 * t220;
	t217 = t212 * t223 + t213 * t221;
	t216 = t212 * t224 - t214 * t220;
	t211 = cos(pkin(9));
	t210 = sin(pkin(9));
	t1 = [t210 * t219 + (-t210 * t225 + t211 * t213) * qJD(1), t216 * t210, 0, 0, 0, 0; -t210 * t218 + (-t210 * t226 - t211 * t215) * qJD(1), -t217 * t210, 0, 0, 0, 0; 0, -t210 * t222, 0, 0, 0, 0; t217, t214 * t224 + t218, 0, 0, 0, 0; t216, -t214 * t223 + t219, 0, 0, 0, 0; 0, -t221, 0, 0, 0, 0; t211 * t219 + (-t210 * t213 - t211 * t225) * qJD(1), t216 * t211, 0, 0, 0, 0; -t211 * t218 + (t210 * t215 - t211 * t226) * qJD(1), -t217 * t211, 0, 0, 0, 0; 0, -t211 * t222, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:11
	% EndTime: 2019-10-10 09:25:11
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (108->40), mult. (385->83), div. (0->0), fcn. (401->8), ass. (0->40)
	t304 = cos(qJ(2));
	t301 = sin(qJ(2));
	t302 = sin(qJ(1));
	t320 = t302 * qJD(2) * t301;
	t305 = cos(qJ(1));
	t324 = qJD(1) * t305;
	t331 = t304 * t324 - t320;
	t298 = sin(pkin(9));
	t299 = cos(pkin(9));
	t300 = sin(qJ(6));
	t303 = cos(qJ(6));
	t316 = t298 * t303 + t299 * t300;
	t330 = qJD(6) * t316;
	t325 = qJD(1) * t302;
	t288 = t331 * t298 - t299 * t325;
	t326 = t305 * t299;
	t293 = t302 * t298 + t304 * t326;
	t289 = t293 * qJD(1) - t299 * t320;
	t328 = t302 * t304;
	t290 = t298 * t328 + t326;
	t327 = t305 * t298;
	t291 = t299 * t328 - t327;
	t329 = t288 * t300 - t289 * t303 + (t290 * t303 + t291 * t300) * qJD(6);
	t323 = qJD(2) * t304;
	t322 = qJD(2) * t305;
	t319 = t301 * t322;
	t315 = t298 * t300 - t299 * t303;
	t314 = t315 * t304;
	t313 = qJD(2) * t316;
	t312 = qJD(2) * t315;
	t309 = qJD(6) * t315;
	t308 = t304 * t312;
	t307 = t304 * t313;
	t306 = -t288 * t303 - t289 * t300 + (t290 * t300 - t291 * t303) * qJD(6);
	t292 = -t302 * t299 + t304 * t327;
	t287 = -t291 * qJD(1) - t299 * t319;
	t286 = -t290 * qJD(1) - t298 * t319;
	t285 = t286 * t303 + t287 * t300 + (-t292 * t300 + t293 * t303) * qJD(6);
	t284 = -t286 * t300 + t287 * t303 + (-t292 * t303 - t293 * t300) * qJD(6);
	t1 = [t306, -t305 * t307 + (t305 * t309 + t316 * t325) * t301, 0, 0, 0, t284; t285, -t302 * t307 + (t302 * t309 - t316 * t324) * t301, 0, 0, 0, -t329; 0, -qJD(6) * t314 - t301 * t313, 0, 0, 0, -qJD(2) * t314 - t301 * t330; t329, t305 * t308 + (t305 * t330 - t315 * t325) * t301, 0, 0, 0, -t285; t284, t302 * t308 + (t302 * t330 + t315 * t324) * t301, 0, 0, 0, t306; 0, t301 * t312 - t304 * t330, 0, 0, 0, t301 * t309 - t307; -t301 * t324 - t302 * t323, -t304 * t325 - t319, 0, 0, 0, 0; -t301 * t325 + t304 * t322, t331, 0, 0, 0, 0; 0, t323, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end