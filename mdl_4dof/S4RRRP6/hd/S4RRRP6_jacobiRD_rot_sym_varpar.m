% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S4RRRP6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRP6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_jacobiRD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:19:24
	% EndTime: 2019-12-31 17:19:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:19:24
	% EndTime: 2019-12-31 17:19:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0; -t31, 0, 0, 0; 0, 0, 0, 0; t31, 0, 0, 0; -t30, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:19:24
	% EndTime: 2019-12-31 17:19:24
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
	t1 = [t32, t29, 0, 0; -t30, -t31, 0, 0; 0, -t39, 0, 0; t31, t30, 0, 0; t29, t32, 0, 0; 0, -t38, 0, 0; -t41, 0, 0, 0; t40, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:19:25
	% EndTime: 2019-12-31 17:19:25
	% DurationCPUTime: 0.07s
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
	t1 = [t228, -t232 * t241 + (t232 * t255 + t234 * t249) * t230, t225, 0; -t226, -t232 * t243 + (t231 * t249 - t232 * t253) * t230, -t227, 0; 0, -t229 * t247 - t246, -t229 * t251 - t245, 0; t227, t236 * t229 + t234 * t245, t226, 0; t225, t237 * t229 + t231 * t245, t228, 0; 0, t229 * t252 - t232 * t247, t229 * t248 - t232 * t251, 0; -t237, -t231 * t254 - t242, 0, 0; t236, t233 * t253 - t244, 0, 0; 0, t251, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:19:25
	% EndTime: 2019-12-31 17:19:25
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (48->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t254 = cos(qJ(3));
	t256 = cos(qJ(1));
	t278 = t254 * t256;
	t253 = sin(qJ(1));
	t277 = qJD(1) * t253;
	t255 = cos(qJ(2));
	t276 = qJD(1) * t255;
	t275 = qJD(1) * t256;
	t252 = sin(qJ(2));
	t274 = qJD(2) * t252;
	t273 = qJD(2) * t255;
	t272 = qJD(2) * t256;
	t251 = sin(qJ(3));
	t271 = qJD(3) * t251;
	t270 = qJD(3) * t252;
	t269 = qJD(3) * t255;
	t268 = t254 * t274;
	t267 = t254 * t270;
	t266 = t253 * t274;
	t265 = t253 * t273;
	t264 = t252 * t272;
	t263 = t255 * t272;
	t262 = -qJD(1) + t269;
	t261 = -qJD(3) + t276;
	t260 = t262 * t251;
	t259 = t252 * t275 + t265;
	t258 = -t252 * t277 + t263;
	t257 = t261 * t253 + t264;
	t250 = -t261 * t278 + (t260 + t268) * t253;
	t249 = t262 * t254 * t253 + (t261 * t256 - t266) * t251;
	t248 = t257 * t254 + t256 * t260;
	t247 = t257 * t251 - t262 * t278;
	t1 = [t250, -t254 * t263 + (t254 * t277 + t256 * t271) * t252, t247, 0; -t248, -t254 * t265 + (t253 * t271 - t254 * t275) * t252, -t249, 0; 0, -t251 * t269 - t268, -t251 * t273 - t267, 0; t249, t258 * t251 + t256 * t267, t248, 0; t247, t259 * t251 + t253 * t267, t250, 0; 0, t251 * t274 - t254 * t269, t251 * t270 - t254 * t273, 0; -t259, -t253 * t276 - t264, 0, 0; t258, t255 * t275 - t266, 0, 0; 0, t273, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,4);
end