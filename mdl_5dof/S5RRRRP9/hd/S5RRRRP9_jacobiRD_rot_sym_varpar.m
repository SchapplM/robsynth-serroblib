% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRP9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRP9_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP9_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:07:14
	% EndTime: 2019-12-31 22:07:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:07:14
	% EndTime: 2019-12-31 22:07:14
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
	% StartTime: 2019-12-31 22:07:14
	% EndTime: 2019-12-31 22:07:14
	% DurationCPUTime: 0.02s
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
	% StartTime: 2019-12-31 22:07:15
	% EndTime: 2019-12-31 22:07:15
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (48->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t232 = cos(qJ(3));
	t234 = cos(qJ(1));
	t256 = t232 * t234;
	t231 = sin(qJ(1));
	t255 = qJD(1) * t231;
	t233 = cos(qJ(2));
	t254 = qJD(1) * t233;
	t253 = qJD(1) * t234;
	t230 = sin(qJ(2));
	t252 = qJD(2) * t230;
	t251 = qJD(2) * t233;
	t250 = qJD(2) * t234;
	t229 = sin(qJ(3));
	t249 = qJD(3) * t229;
	t248 = qJD(3) * t230;
	t247 = qJD(3) * t233;
	t246 = t232 * t252;
	t245 = t232 * t248;
	t244 = t231 * t252;
	t243 = t231 * t251;
	t242 = t230 * t250;
	t241 = t233 * t250;
	t240 = -qJD(1) + t247;
	t239 = -qJD(3) + t254;
	t238 = t240 * t229;
	t237 = t230 * t253 + t243;
	t236 = -t230 * t255 + t241;
	t235 = t239 * t231 + t242;
	t228 = -t239 * t256 + (t238 + t246) * t231;
	t227 = t240 * t232 * t231 + (t239 * t234 - t244) * t229;
	t226 = t235 * t232 + t234 * t238;
	t225 = t235 * t229 - t240 * t256;
	t1 = [t228, -t232 * t241 + (t232 * t255 + t234 * t249) * t230, t225, 0, 0; -t226, -t232 * t243 + (t231 * t249 - t232 * t253) * t230, -t227, 0, 0; 0, -t229 * t247 - t246, -t229 * t251 - t245, 0, 0; t227, t236 * t229 + t234 * t245, t226, 0, 0; t225, t237 * t229 + t231 * t245, t228, 0, 0; 0, t229 * t252 - t232 * t247, t229 * t248 - t232 * t251, 0, 0; -t237, -t231 * t254 - t242, 0, 0, 0; t236, t233 * t253 - t244, 0, 0, 0; 0, t251, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:07:15
	% EndTime: 2019-12-31 22:07:15
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (184->30), mult. (233->48), div. (0->0), fcn. (233->6), ass. (0->35)
	t288 = qJD(3) + qJD(4);
	t290 = sin(qJ(2));
	t313 = t288 * t290;
	t292 = cos(qJ(2));
	t312 = t288 * t292;
	t291 = sin(qJ(1));
	t311 = qJD(1) * t291;
	t310 = qJD(1) * t292;
	t293 = cos(qJ(1));
	t309 = qJD(1) * t293;
	t308 = qJD(2) * t290;
	t307 = qJD(2) * t292;
	t306 = qJD(2) * t293;
	t289 = qJ(3) + qJ(4);
	t286 = sin(t289);
	t305 = t286 * t313;
	t287 = cos(t289);
	t304 = t287 * t313;
	t303 = t291 * t288 * t287;
	t302 = t293 * t288 * t286;
	t301 = t291 * t308;
	t300 = t290 * t306;
	t299 = -qJD(1) + t312;
	t298 = -t288 + t310;
	t297 = t286 * t299;
	t296 = t290 * t309 + t291 * t307;
	t295 = t290 * t311 - t292 * t306;
	t294 = t298 * t291 + t300;
	t282 = -t287 * t307 + t305;
	t281 = -t286 * t307 - t304;
	t280 = t291 * t297 + (-t298 * t293 + t301) * t287;
	t279 = -t286 * t301 - t302 - t287 * t311 + (t286 * t309 + t303) * t292;
	t278 = t294 * t287 + t293 * t297;
	t277 = -t299 * t293 * t287 + t294 * t286;
	t1 = [t280, t295 * t287 + t290 * t302, t277, t277, 0; -t278, -t296 * t287 + t291 * t305, -t279, -t279, 0; 0, -t286 * t312 - t287 * t308, t281, t281, 0; t279, -t295 * t286 + t293 * t304, t278, t278, 0; t277, t296 * t286 + t290 * t303, t280, t280, 0; 0, t286 * t308 - t287 * t312, t282, t282, 0; -t296, -t291 * t310 - t300, 0, 0, 0; -t295, t292 * t309 - t301, 0, 0, 0; 0, t307, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:07:16
	% EndTime: 2019-12-31 22:07:16
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (184->25), mult. (233->44), div. (0->0), fcn. (233->6), ass. (0->34)
	t365 = sin(qJ(1));
	t362 = qJD(3) + qJD(4);
	t366 = cos(qJ(2));
	t382 = qJD(1) * t366;
	t372 = -t362 + t382;
	t364 = sin(qJ(2));
	t367 = cos(qJ(1));
	t378 = qJD(2) * t367;
	t374 = t364 * t378;
	t386 = t372 * t365 + t374;
	t380 = qJD(2) * t364;
	t375 = t365 * t380;
	t385 = t372 * t367 - t375;
	t384 = t362 * t364;
	t383 = t362 * t366;
	t381 = qJD(1) * t367;
	t379 = qJD(2) * t366;
	t363 = qJ(3) + qJ(4);
	t360 = sin(t363);
	t377 = t360 * t384;
	t361 = cos(t363);
	t376 = t361 * t384;
	t373 = qJD(1) - t383;
	t371 = t373 * t365;
	t370 = t373 * t367;
	t369 = -t364 * t381 - t365 * t379;
	t368 = qJD(1) * t365 * t364 - t366 * t378;
	t359 = t361 * t379 - t377;
	t358 = -t360 * t379 - t376;
	t357 = t360 * t371 + t385 * t361;
	t356 = -t385 * t360 + t361 * t371;
	t355 = t360 * t370 - t386 * t361;
	t354 = t386 * t360 + t361 * t370;
	t1 = [-t357, t368 * t361 + t367 * t377, t354, t354, 0; t355, t369 * t361 + t365 * t377, t356, t356, 0; 0, -t360 * t383 - t361 * t380, t358, t358, 0; t369, -t365 * t382 - t374, 0, 0, 0; -t368, t366 * t381 - t375, 0, 0, 0; 0, t379, 0, 0, 0; t356, t368 * t360 - t367 * t376, t355, t355, 0; -t354, t369 * t360 - t365 * t376, t357, t357, 0; 0, -t360 * t380 + t361 * t383, t359, t359, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end