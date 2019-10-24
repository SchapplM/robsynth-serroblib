% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPPRR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:17
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PPPRR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->8), mult. (42->25), div. (0->0), fcn. (46->8), ass. (0->14)
	t84 = sin(pkin(8));
	t89 = sin(qJ(4));
	t93 = t84 * t89;
	t90 = cos(qJ(4));
	t92 = t84 * t90;
	t86 = cos(pkin(9));
	t87 = cos(pkin(8));
	t91 = t86 * t87;
	t88 = cos(pkin(7));
	t85 = sin(pkin(7));
	t83 = sin(pkin(9));
	t82 = t85 * t83 + t88 * t91;
	t81 = -t88 * t83 + t85 * t91;
	t1 = [0, 0, 0, (-t82 * t90 - t88 * t93) * qJD(4), 0; 0, 0, 0, (-t81 * t90 - t85 * t93) * qJD(4), 0; 0, 0, 0, (-t86 * t92 + t87 * t89) * qJD(4), 0; 0, 0, 0, (t82 * t89 - t88 * t92) * qJD(4), 0; 0, 0, 0, (t81 * t89 - t85 * t92) * qJD(4), 0; 0, 0, 0, (t86 * t93 + t87 * t90) * qJD(4), 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:42
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (67->28), mult. (227->68), div. (0->0), fcn. (261->10), ass. (0->35)
	t265 = sin(pkin(9));
	t266 = sin(pkin(8));
	t282 = t265 * t266;
	t272 = sin(qJ(4));
	t281 = t266 * t272;
	t274 = cos(qJ(4));
	t280 = t266 * t274;
	t267 = sin(pkin(7));
	t269 = cos(pkin(8));
	t279 = t267 * t269;
	t270 = cos(pkin(7));
	t278 = t270 * t265;
	t268 = cos(pkin(9));
	t277 = t270 * t268;
	t271 = sin(qJ(5));
	t276 = qJD(5) * t271;
	t273 = cos(qJ(5));
	t275 = qJD(5) * t273;
	t260 = t268 * t279 - t278;
	t253 = -t260 * t272 + t267 * t280;
	t254 = t260 * t274 + t267 * t281;
	t262 = t265 * t267 + t269 * t277;
	t255 = -t262 * t272 + t270 * t280;
	t256 = t262 * t274 + t270 * t281;
	t264 = t268 * t280 - t269 * t272;
	t263 = -t268 * t281 - t269 * t274;
	t261 = -t267 * t268 + t269 * t278;
	t259 = t265 * t279 + t277;
	t258 = t264 * qJD(4);
	t257 = t263 * qJD(4);
	t252 = t256 * qJD(4);
	t251 = t255 * qJD(4);
	t250 = t254 * qJD(4);
	t249 = t253 * qJD(4);
	t1 = [0, 0, 0, -t252 * t273 - t255 * t276, -t251 * t271 + (-t256 * t273 - t261 * t271) * qJD(5); 0, 0, 0, -t250 * t273 - t253 * t276, -t249 * t271 + (-t254 * t273 - t259 * t271) * qJD(5); 0, 0, 0, -t258 * t273 - t263 * t276, -t257 * t271 + (-t264 * t273 - t271 * t282) * qJD(5); 0, 0, 0, t252 * t271 - t255 * t275, -t251 * t273 + (t256 * t271 - t261 * t273) * qJD(5); 0, 0, 0, t250 * t271 - t253 * t275, -t249 * t273 + (t254 * t271 - t259 * t273) * qJD(5); 0, 0, 0, t258 * t271 - t263 * t275, -t257 * t273 + (t264 * t271 - t273 * t282) * qJD(5); 0, 0, 0, t251, 0; 0, 0, 0, t249, 0; 0, 0, 0, t257, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end