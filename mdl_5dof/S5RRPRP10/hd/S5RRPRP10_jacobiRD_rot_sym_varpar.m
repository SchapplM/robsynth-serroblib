% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRP10
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRP10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:18
	% EndTime: 2019-12-29 18:59:18
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:18
	% EndTime: 2019-12-29 18:59:18
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:13
	% EndTime: 2019-12-29 18:59:13
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-12-29 18:59:19
	% EndTime: 2019-12-29 18:59:19
	% DurationCPUTime: 0.06s
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
	t1 = [-t168, 0, 0, 0, 0; t167, 0, 0, 0, 0; 0, 0, 0, 0, 0; t159, t156, 0, 0, 0; t157, t158, 0, 0, 0; 0, t166, 0, 0, 0; -t158, -t157, 0, 0, 0; t156, t159, 0, 0, 0; 0, t165, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:15
	% EndTime: 2019-12-29 18:59:15
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (49->25), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->33)
	t237 = cos(qJ(4));
	t235 = sin(qJ(2));
	t254 = qJD(4) * t235;
	t246 = qJD(1) + t254;
	t261 = t237 * t246;
	t234 = sin(qJ(4));
	t260 = t246 * t234;
	t236 = sin(qJ(1));
	t259 = qJD(1) * t236;
	t239 = cos(qJ(1));
	t258 = qJD(1) * t239;
	t257 = qJD(2) * t235;
	t238 = cos(qJ(2));
	t256 = qJD(2) * t238;
	t255 = qJD(2) * t239;
	t253 = qJD(4) * t238;
	t252 = t238 * t258;
	t251 = t237 * t256;
	t250 = t236 * t256;
	t249 = t234 * t253;
	t248 = t237 * t253;
	t247 = t238 * t255;
	t245 = -qJD(1) * t235 - qJD(4);
	t244 = t245 * t239;
	t243 = -t236 * t257 + t252;
	t242 = -t235 * t255 - t238 * t259;
	t241 = t237 * t257 + t249;
	t240 = t245 * t236 + t247;
	t233 = t240 * t234 + t239 * t261;
	t232 = t240 * t237 - t239 * t260;
	t231 = -t236 * t261 + (t244 - t250) * t234;
	t230 = t237 * t244 + (-t251 + t260) * t236;
	t1 = [t231, t242 * t234 + t239 * t248, 0, t232, 0; t233, t243 * t234 + t236 * t248, 0, -t230, 0; 0, t234 * t256 + t237 * t254, 0, t241, 0; t230, t242 * t237 - t239 * t249, 0, -t233, 0; t232, -t241 * t236 + t237 * t252, 0, t231, 0; 0, -t234 * t254 + t251, 0, -t234 * t257 + t248, 0; -t243, t235 * t259 - t247, 0, 0, 0; t242, -t235 * t258 - t250, 0, 0, 0; 0, -t257, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:15
	% EndTime: 2019-12-29 18:59:15
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (49->25), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->33)
	t253 = cos(qJ(4));
	t251 = sin(qJ(2));
	t270 = qJD(4) * t251;
	t262 = qJD(1) + t270;
	t277 = t253 * t262;
	t250 = sin(qJ(4));
	t276 = t262 * t250;
	t252 = sin(qJ(1));
	t275 = qJD(1) * t252;
	t255 = cos(qJ(1));
	t274 = qJD(1) * t255;
	t273 = qJD(2) * t251;
	t254 = cos(qJ(2));
	t272 = qJD(2) * t254;
	t271 = qJD(2) * t255;
	t269 = qJD(4) * t254;
	t268 = t254 * t274;
	t267 = t253 * t272;
	t266 = t252 * t272;
	t265 = t250 * t269;
	t264 = t253 * t269;
	t263 = t254 * t271;
	t261 = -qJD(1) * t251 - qJD(4);
	t260 = t261 * t255;
	t259 = -t252 * t273 + t268;
	t258 = -t251 * t271 - t254 * t275;
	t257 = t253 * t273 + t265;
	t256 = t261 * t252 + t263;
	t249 = t256 * t250 + t255 * t277;
	t248 = t256 * t253 - t255 * t276;
	t247 = -t252 * t277 + (t260 - t266) * t250;
	t246 = t253 * t260 + (-t267 + t276) * t252;
	t1 = [t247, t258 * t250 + t255 * t264, 0, t248, 0; t249, t259 * t250 + t252 * t264, 0, -t246, 0; 0, t250 * t272 + t253 * t270, 0, t257, 0; t246, t258 * t253 - t255 * t265, 0, -t249, 0; t248, -t257 * t252 + t253 * t268, 0, t247, 0; 0, -t250 * t270 + t267, 0, -t250 * t273 + t264, 0; -t259, t251 * t275 - t263, 0, 0, 0; t258, -t251 * t274 - t266, 0, 0, 0; 0, -t273, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end