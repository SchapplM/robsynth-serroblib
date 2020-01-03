% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPPR10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:31
	% EndTime: 2019-12-31 19:45:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:32
	% EndTime: 2019-12-31 19:45:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:32
	% EndTime: 2019-12-31 19:45:32
	% DurationCPUTime: 0.03s
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
	t1 = [t32, t29, 0, 0, 0; -t30, -t31, 0, 0, 0; 0, -t39, 0, 0, 0; t31, t30, 0, 0, 0; t29, t32, 0, 0, 0; 0, -t38, 0, 0, 0; -t41, 0, 0, 0, 0; t40, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:32
	% EndTime: 2019-12-31 19:45:33
	% DurationCPUTime: 0.05s
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
	t178 = cos(pkin(8));
	t177 = sin(pkin(8));
	t1 = [t178 * t186 + (-t177 * t180 - t178 * t192) * qJD(1), t183 * t178, 0, 0, 0; -t178 * t185 + (t177 * t182 - t178 * t193) * qJD(1), -t184 * t178, 0, 0, 0; 0, -t178 * t189, 0, 0, 0; -t177 * t186 + (t177 * t192 - t178 * t180) * qJD(1), -t183 * t177, 0, 0, 0; t177 * t185 + (t177 * t193 + t178 * t182) * qJD(1), t184 * t177, 0, 0, 0; 0, t177 * t189, 0, 0, 0; -t184, -t181 * t191 - t185, 0, 0, 0; -t183, t181 * t190 - t186, 0, 0, 0; 0, t188, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:33
	% EndTime: 2019-12-31 19:45:33
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (18->15), mult. (77->37), div. (0->0), fcn. (77->6), ass. (0->18)
	t215 = sin(qJ(1));
	t216 = cos(qJ(2));
	t228 = t215 * t216;
	t217 = cos(qJ(1));
	t227 = t216 * t217;
	t226 = qJD(1) * t215;
	t225 = qJD(1) * t217;
	t214 = sin(qJ(2));
	t224 = qJD(2) * t214;
	t223 = qJD(2) * t216;
	t222 = qJD(2) * t217;
	t221 = t215 * t224;
	t220 = t214 * t222;
	t219 = -t214 * t225 - t215 * t223;
	t218 = t214 * t226 - t216 * t222;
	t213 = cos(pkin(8));
	t212 = sin(pkin(8));
	t1 = [t213 * t221 + (-t212 * t215 - t213 * t227) * qJD(1), t218 * t213, 0, 0, 0; -t213 * t220 + (t212 * t217 - t213 * t228) * qJD(1), t219 * t213, 0, 0, 0; 0, -t213 * t224, 0, 0, 0; t219, -t216 * t226 - t220, 0, 0, 0; -t218, t216 * t225 - t221, 0, 0, 0; 0, t223, 0, 0, 0; t212 * t221 + (-t212 * t227 + t213 * t215) * qJD(1), t218 * t212, 0, 0, 0; -t212 * t220 + (-t212 * t228 - t213 * t217) * qJD(1), t219 * t212, 0, 0, 0; 0, -t212 * t224, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:33
	% EndTime: 2019-12-31 19:45:33
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (109->41), mult. (385->83), div. (0->0), fcn. (401->8), ass. (0->40)
	t302 = cos(qJ(2));
	t299 = sin(qJ(2));
	t300 = sin(qJ(1));
	t318 = t300 * qJD(2) * t299;
	t303 = cos(qJ(1));
	t322 = qJD(1) * t303;
	t329 = -t302 * t322 + t318;
	t296 = sin(pkin(8));
	t297 = cos(pkin(8));
	t298 = sin(qJ(5));
	t301 = cos(qJ(5));
	t313 = t296 * t298 + t297 * t301;
	t328 = qJD(5) * t313;
	t323 = qJD(1) * t300;
	t286 = -t329 * t296 - t297 * t323;
	t324 = t303 * t297;
	t291 = t300 * t296 + t302 * t324;
	t287 = t291 * qJD(1) - t297 * t318;
	t326 = t300 * t302;
	t288 = t296 * t326 + t324;
	t325 = t303 * t296;
	t289 = t297 * t326 - t325;
	t327 = -t286 * t301 + t287 * t298 + (t288 * t298 + t289 * t301) * qJD(5);
	t321 = qJD(2) * t302;
	t320 = qJD(2) * t303;
	t317 = t299 * t320;
	t314 = t296 * t301 - t297 * t298;
	t312 = t314 * t302;
	t311 = qJD(2) * t314;
	t310 = qJD(2) * t313;
	t309 = qJD(5) * t314;
	t306 = t302 * t311;
	t305 = t302 * t310;
	t304 = -t286 * t298 - t287 * t301 + (-t288 * t301 + t289 * t298) * qJD(5);
	t290 = -t300 * t297 + t302 * t325;
	t285 = -t289 * qJD(1) - t297 * t317;
	t284 = -t288 * qJD(1) - t296 * t317;
	t283 = t284 * t298 + t285 * t301 + (t290 * t301 - t291 * t298) * qJD(5);
	t282 = t284 * t301 - t285 * t298 + (-t290 * t298 - t291 * t301) * qJD(5);
	t1 = [t304, -t303 * t305 + (-t303 * t309 + t313 * t323) * t299, 0, 0, t282; t283, -t300 * t305 + (-t300 * t309 - t313 * t322) * t299, 0, 0, -t327; 0, qJD(5) * t312 - t299 * t310, 0, 0, qJD(2) * t312 - t299 * t328; t327, -t303 * t306 + (t303 * t328 + t314 * t323) * t299, 0, 0, -t283; t282, -t300 * t306 + (t300 * t328 - t314 * t322) * t299, 0, 0, t304; 0, -t299 * t311 - t302 * t328, 0, 0, -t299 * t309 - t305; t299 * t322 + t300 * t321, t302 * t323 + t317, 0, 0, 0; t299 * t323 - t302 * t320, t329, 0, 0, 0; 0, -t321, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end