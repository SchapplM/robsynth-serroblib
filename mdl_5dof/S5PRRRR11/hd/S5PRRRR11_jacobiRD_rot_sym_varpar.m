% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRR11
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
%   Siehe auch: S5PRRRR11_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRRRR11_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR11_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	% Symbolic code from jacobiRD_rot_1_floatb_twist_matlab.m not found
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = pkin(10) + qJ(2);
	t37 = qJD(2) * sin(t35);
	t36 = qJD(2) * cos(t35);
	t1 = [0, -t36, 0, 0, 0; 0, -t37, 0, 0, 0; 0, 0, 0, 0, 0; 0, t37, 0, 0, 0; 0, -t36, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:37
	% EndTime: 2024-09-27 21:46:37
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (61->14), mult. (88->24), div. (0->0), fcn. (88->6), ass. (0->20)
	t161 = cos(pkin(5));
	t162 = sin(qJ(3));
	t171 = t161 * t162;
	t163 = cos(qJ(3));
	t170 = t161 * t163;
	t160 = sin(pkin(5));
	t169 = qJD(2) * t160;
	t168 = qJD(3) * t160;
	t159 = pkin(10) + qJ(2);
	t157 = sin(t159);
	t158 = cos(t159);
	t167 = t157 * t162 - t158 * t170;
	t166 = t157 * t163 + t158 * t171;
	t165 = t157 * t170 + t158 * t162;
	t164 = t157 * t171 - t158 * t163;
	t156 = t164 * qJD(2) + t167 * qJD(3);
	t155 = t165 * qJD(2) + t166 * qJD(3);
	t154 = t166 * qJD(2) + t165 * qJD(3);
	t153 = t167 * qJD(2) + t164 * qJD(3);
	t1 = [0, t156, t153, 0, 0; 0, -t154, -t155, 0, 0; 0, 0, -t162 * t168, 0, 0; 0, t155, t154, 0, 0; 0, t153, t156, 0, 0; 0, 0, -t163 * t168, 0, 0; 0, -t157 * t169, 0, 0, 0; 0, t158 * t169, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:37
	% EndTime: 2024-09-27 21:46:37
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (295->26), mult. (176->37), div. (0->0), fcn. (132->9), ass. (0->30)
	t278 = qJ(3) + qJ(4);
	t273 = pkin(5) - t278;
	t288 = cos(t273);
	t272 = pkin(5) + t278;
	t290 = cos(t272) / 0.2e1;
	t261 = t288 / 0.2e1 + t290;
	t277 = qJD(3) + qJD(4);
	t259 = t261 * t277;
	t268 = sin(t273);
	t287 = sin(t272);
	t260 = t287 / 0.2e1 - t268 / 0.2e1;
	t276 = pkin(10) + qJ(2);
	t270 = sin(t276);
	t271 = cos(t276);
	t275 = cos(t278);
	t274 = sin(t278);
	t286 = t274 * t277;
	t253 = (t260 * t271 + t270 * t275) * qJD(2) + t270 * t259 + t271 * t286;
	t289 = -t277 / 0.2e1;
	t285 = t277 * t268;
	t284 = t277 * t275;
	t283 = qJD(2) * sin(pkin(5));
	t263 = t287 * t289;
	t282 = t285 / 0.2e1 + t263;
	t255 = -t271 * t259 + t270 * t286 + (t260 * t270 - t271 * t275) * qJD(2);
	t258 = t277 * t290 + t288 * t289;
	t257 = t263 - t285 / 0.2e1;
	t254 = t270 * t284 - t271 * t282 + (t261 * t270 + t271 * t274) * qJD(2);
	t252 = -t271 * t284 - t270 * t282 + (-t271 * t261 + t270 * t274) * qJD(2);
	t1 = [0, t255, t252, t252, 0; 0, -t253, -t254, -t254, 0; 0, 0, t258, t258, 0; 0, t254, t253, t253, 0; 0, t252, t255, t255, 0; 0, 0, t257, t257, 0; 0, -t270 * t283, 0, 0, 0; 0, t271 * t283, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:37
	% EndTime: 2024-09-27 21:46:37
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (565->26), mult. (236->34), div. (0->0), fcn. (176->9), ass. (0->28)
	t308 = qJ(3) + qJ(4) + qJ(5);
	t302 = pkin(5) - t308;
	t321 = cos(t302);
	t301 = pkin(5) + t308;
	t322 = cos(t301) / 0.2e1;
	t294 = t321 / 0.2e1 + t322;
	t307 = qJD(3) + qJD(4) + qJD(5);
	t292 = t294 * t307;
	t304 = cos(t308);
	t309 = pkin(10) + qJ(2);
	t305 = sin(t309);
	t306 = cos(t309);
	t315 = -sin(t301) / 0.2e1;
	t320 = sin(t302);
	t313 = t320 / 0.2e1 + t315;
	t303 = sin(t308);
	t319 = t303 * t306;
	t286 = (t304 * t305 - t306 * t313) * qJD(2) + t305 * t292 + t307 * t319;
	t318 = t305 * t303;
	t317 = t307 * t304;
	t316 = qJD(2) * sin(pkin(5));
	t311 = t313 * t307;
	t288 = -t306 * t292 + t307 * t318 + (-t306 * t304 - t305 * t313) * qJD(2);
	t291 = (t322 - t321 / 0.2e1) * t307;
	t290 = (t315 - t320 / 0.2e1) * t307;
	t287 = t305 * t317 - t306 * t311 + (t294 * t305 + t319) * qJD(2);
	t285 = -t305 * t311 - t306 * t317 + (-t306 * t294 + t318) * qJD(2);
	t1 = [0, t288, t285, t285, t285; 0, -t286, -t287, -t287, -t287; 0, 0, t291, t291, t291; 0, t287, t286, t286, t286; 0, t285, t288, t288, t288; 0, 0, t290, t290, t290; 0, -t305 * t316, 0, 0, 0; 0, t306 * t316, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end