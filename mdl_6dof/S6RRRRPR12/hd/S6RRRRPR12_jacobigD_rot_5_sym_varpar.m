% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR12
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR12_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:59
% EndTime: 2019-02-26 22:36:59
% DurationCPUTime: 0.13s
% Computational Cost: add. (38->24), mult. (136->54), div. (0->0), fcn. (140->10), ass. (0->31)
t224 = sin(pkin(7));
t225 = sin(pkin(6));
t251 = t225 * t224;
t226 = cos(pkin(7));
t231 = cos(qJ(3));
t250 = t226 * t231;
t228 = sin(qJ(3));
t232 = cos(qJ(2));
t249 = t228 * t232;
t229 = sin(qJ(2));
t248 = t229 * t231;
t230 = sin(qJ(1));
t247 = t230 * t229;
t246 = t230 * t232;
t233 = cos(qJ(1));
t245 = t233 * t229;
t244 = t233 * t232;
t243 = qJD(1) * t225;
t242 = qJD(2) * t228;
t241 = t230 * t251;
t240 = t233 * t251;
t239 = t230 * t243;
t238 = t233 * t243;
t227 = cos(pkin(6));
t237 = t227 * t244 - t247;
t236 = -t227 * t246 - t245;
t235 = t227 * t245 + t246;
t234 = t227 * t247 - t244;
t223 = qJD(1) * t236 - qJD(2) * t235;
t222 = -qJD(1) * t237 + qJD(2) * t234;
t1 = [0, t238, -t222 * t224 + t226 * t238, -t222 * t250 + t236 * t242 + (-t228 * t235 - t231 * t240) * qJD(1) + (-t234 * t231 + (t226 * t236 + t241) * t228) * qJD(3), 0, 0; 0, t239, -t223 * t224 + t226 * t239, -t223 * t250 + t237 * t242 + (-t228 * t234 - t231 * t241) * qJD(1) + (t235 * t231 + (t226 * t237 - t240) * t228) * qJD(3), 0, 0; 0, 0, qJD(2) * t229 * t251, t227 * t224 * qJD(3) * t228 + ((t226 * t249 + t248) * qJD(3) + (t226 * t248 + t249) * qJD(2)) * t225, 0, 0;];
JgD_rot  = t1;
