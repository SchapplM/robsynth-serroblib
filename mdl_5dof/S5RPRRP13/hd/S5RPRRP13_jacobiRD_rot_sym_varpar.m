% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRP13
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRP13_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP13_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:00:06
	% EndTime: 2019-12-31 19:00:06
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:00:06
	% EndTime: 2019-12-31 19:00:07
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
	% StartTime: 2019-12-31 19:00:06
	% EndTime: 2019-12-31 19:00:06
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t11, 0, 0, 0, 0; t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0; t11, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:00:06
	% EndTime: 2019-12-31 19:00:06
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t35 = sin(qJ(1));
	t42 = qJD(1) * t35;
	t37 = cos(qJ(1));
	t41 = qJD(1) * t37;
	t34 = sin(qJ(3));
	t40 = qJD(3) * t34;
	t36 = cos(qJ(3));
	t39 = qJD(3) * t36;
	t38 = qJD(3) * t37;
	t33 = -t35 * t40 + t36 * t41;
	t32 = t34 * t41 + t35 * t39;
	t31 = t34 * t38 + t36 * t42;
	t30 = -t34 * t42 + t36 * t38;
	t1 = [t30, 0, t33, 0, 0; t32, 0, t31, 0, 0; 0, 0, -t39, 0, 0; -t31, 0, -t32, 0, 0; t33, 0, t30, 0, 0; 0, 0, t40, 0, 0; -t41, 0, 0, 0, 0; -t42, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:00:07
	% EndTime: 2019-12-31 19:00:08
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (49->24), mult. (173->45), div. (0->0), fcn. (173->6), ass. (0->32)
	t239 = cos(qJ(1));
	t235 = sin(qJ(3));
	t244 = qJD(1) * t235 + qJD(4);
	t260 = t244 * t239;
	t236 = sin(qJ(1));
	t238 = cos(qJ(3));
	t254 = qJD(3) * t239;
	t246 = t238 * t254;
	t259 = t244 * t236 - t246;
	t258 = qJD(1) * t236;
	t257 = qJD(1) * t239;
	t256 = qJD(3) * t235;
	t255 = qJD(3) * t238;
	t253 = qJD(4) * t235;
	t252 = qJD(4) * t238;
	t251 = t238 * t257;
	t237 = cos(qJ(4));
	t250 = t237 * t255;
	t249 = t236 * t255;
	t234 = sin(qJ(4));
	t248 = t234 * t252;
	t247 = t237 * t252;
	t245 = -qJD(1) - t253;
	t243 = t245 * t239;
	t242 = t236 * t256 - t251;
	t241 = t235 * t254 + t238 * t258;
	t240 = t237 * t256 + t248;
	t233 = t237 * t260 + (t245 * t234 + t250) * t236;
	t232 = t245 * t237 * t236 + (-t249 - t260) * t234;
	t231 = t234 * t243 - t259 * t237;
	t230 = t259 * t234 + t237 * t243;
	t1 = [t231, 0, -t240 * t236 + t237 * t251, t232, 0; t233, 0, t241 * t237 + t239 * t248, -t230, 0; 0, 0, t234 * t253 - t250, t234 * t256 - t247, 0; t230, 0, t242 * t234 - t236 * t247, -t233, 0; t232, 0, -t241 * t234 + t239 * t247, t231, 0; 0, 0, t234 * t255 + t237 * t253, t240, 0; t241, 0, t235 * t257 + t249, 0, 0; t242, 0, t235 * t258 - t246, 0, 0; 0, 0, -t256, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:00:08
	% EndTime: 2019-12-31 19:00:08
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (49->25), mult. (173->45), div. (0->0), fcn. (173->6), ass. (0->32)
	t281 = sin(qJ(1));
	t280 = sin(qJ(3));
	t290 = qJD(1) * t280 + qJD(4);
	t283 = cos(qJ(3));
	t284 = cos(qJ(1));
	t300 = qJD(3) * t284;
	t292 = t283 * t300;
	t305 = t290 * t281 - t292;
	t304 = qJD(1) * t281;
	t303 = qJD(1) * t284;
	t302 = qJD(3) * t280;
	t301 = qJD(3) * t283;
	t299 = qJD(4) * t280;
	t298 = qJD(4) * t283;
	t297 = t283 * t303;
	t282 = cos(qJ(4));
	t296 = t282 * t301;
	t295 = t281 * t301;
	t279 = sin(qJ(4));
	t294 = t279 * t298;
	t293 = t282 * t298;
	t291 = qJD(1) + t299;
	t289 = t291 * t284;
	t288 = t290 * t284;
	t287 = -t281 * t302 + t297;
	t286 = t280 * t300 + t283 * t304;
	t285 = -t282 * t302 - t294;
	t278 = t282 * t288 + (-t291 * t279 + t296) * t281;
	t277 = t291 * t282 * t281 + (t288 + t295) * t279;
	t276 = t279 * t289 + t305 * t282;
	t275 = -t305 * t279 + t282 * t289;
	t1 = [-t276, 0, t285 * t281 + t282 * t297, -t277, 0; t278, 0, t286 * t282 + t284 * t294, t275, 0; 0, 0, t279 * t299 - t296, t279 * t302 - t293, 0; t286, 0, t280 * t303 + t295, 0, 0; -t287, 0, t280 * t304 - t292, 0, 0; 0, 0, -t302, 0, 0; t275, 0, t287 * t279 + t281 * t293, t278, 0; t277, 0, t286 * t279 - t284 * t293, t276, 0; 0, 0, -t279 * t301 - t282 * t299, t285, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end