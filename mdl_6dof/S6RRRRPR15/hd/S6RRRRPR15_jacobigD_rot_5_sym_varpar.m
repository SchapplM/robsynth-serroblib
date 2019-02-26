% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR15_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobigD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:58
% EndTime: 2019-02-26 22:38:58
% DurationCPUTime: 0.08s
% Computational Cost: add. (38->24), mult. (136->54), div. (0->0), fcn. (140->10), ass. (0->31)
t240 = sin(pkin(7));
t241 = sin(pkin(6));
t267 = t241 * t240;
t242 = cos(pkin(7));
t247 = cos(qJ(3));
t266 = t242 * t247;
t244 = sin(qJ(3));
t248 = cos(qJ(2));
t265 = t244 * t248;
t245 = sin(qJ(2));
t264 = t245 * t247;
t246 = sin(qJ(1));
t263 = t246 * t245;
t262 = t246 * t248;
t249 = cos(qJ(1));
t261 = t249 * t245;
t260 = t249 * t248;
t259 = qJD(1) * t241;
t258 = qJD(2) * t244;
t257 = t246 * t267;
t256 = t249 * t267;
t255 = t246 * t259;
t254 = t249 * t259;
t243 = cos(pkin(6));
t253 = t243 * t260 - t263;
t252 = -t243 * t262 - t261;
t251 = t243 * t261 + t262;
t250 = t243 * t263 - t260;
t239 = t252 * qJD(1) - t251 * qJD(2);
t238 = -t253 * qJD(1) + t250 * qJD(2);
t1 = [0, t254, -t238 * t240 + t242 * t254, -t238 * t266 + t252 * t258 + (-t251 * t244 - t247 * t256) * qJD(1) + (-t250 * t247 + (t252 * t242 + t257) * t244) * qJD(3), 0, 0; 0, t255, -t239 * t240 + t242 * t255, -t239 * t266 + t253 * t258 + (-t250 * t244 - t247 * t257) * qJD(1) + (t251 * t247 + (t253 * t242 - t256) * t244) * qJD(3), 0, 0; 0, 0, qJD(2) * t245 * t267, t243 * t240 * qJD(3) * t244 + ((t242 * t265 + t264) * qJD(3) + (t242 * t264 + t265) * qJD(2)) * t241, 0, 0;];
JgD_rot  = t1;
