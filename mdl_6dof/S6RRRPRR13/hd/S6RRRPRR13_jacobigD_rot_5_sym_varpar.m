% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR13_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:05
% EndTime: 2019-02-26 22:23:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (38->24), mult. (136->54), div. (0->0), fcn. (140->10), ass. (0->31)
t223 = sin(pkin(7));
t224 = sin(pkin(6));
t250 = t224 * t223;
t225 = cos(pkin(7));
t230 = cos(qJ(3));
t249 = t225 * t230;
t227 = sin(qJ(3));
t231 = cos(qJ(2));
t248 = t227 * t231;
t228 = sin(qJ(2));
t247 = t228 * t230;
t229 = sin(qJ(1));
t246 = t229 * t228;
t245 = t229 * t231;
t232 = cos(qJ(1));
t244 = t232 * t228;
t243 = t232 * t231;
t242 = qJD(1) * t224;
t241 = qJD(2) * t227;
t240 = t229 * t250;
t239 = t232 * t250;
t238 = t229 * t242;
t237 = t232 * t242;
t226 = cos(pkin(6));
t236 = t226 * t243 - t246;
t235 = -t226 * t245 - t244;
t234 = t226 * t244 + t245;
t233 = t226 * t246 - t243;
t222 = t235 * qJD(1) - t234 * qJD(2);
t221 = -t236 * qJD(1) + t233 * qJD(2);
t1 = [0, t237, -t221 * t223 + t225 * t237, 0, -t221 * t249 + t235 * t241 + (-t234 * t227 - t230 * t239) * qJD(1) + (-t233 * t230 + (t235 * t225 + t240) * t227) * qJD(3), 0; 0, t238, -t222 * t223 + t225 * t238, 0, -t222 * t249 + t236 * t241 + (-t233 * t227 - t230 * t240) * qJD(1) + (t234 * t230 + (t236 * t225 - t239) * t227) * qJD(3), 0; 0, 0, qJD(2) * t228 * t250, 0, t226 * t223 * qJD(3) * t227 + ((t225 * t248 + t247) * qJD(3) + (t225 * t247 + t248) * qJD(2)) * t224, 0;];
JgD_rot  = t1;
