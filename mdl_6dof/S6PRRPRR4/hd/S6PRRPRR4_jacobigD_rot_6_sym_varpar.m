% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR4_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:53
% EndTime: 2019-02-26 20:05:53
% DurationCPUTime: 0.10s
% Computational Cost: add. (44->21), mult. (144->44), div. (0->0), fcn. (158->10), ass. (0->23)
t250 = qJD(5) - qJD(3);
t225 = sin(pkin(6));
t229 = sin(qJ(3));
t246 = t225 * t229;
t230 = sin(qJ(2));
t245 = t225 * t230;
t232 = cos(qJ(3));
t244 = t225 * t232;
t227 = cos(pkin(6));
t243 = t227 * t230;
t233 = cos(qJ(2));
t242 = t227 * t233;
t241 = qJD(2) * t245;
t224 = sin(pkin(11));
t226 = cos(pkin(11));
t222 = t224 * t233 + t226 * t243;
t223 = -t224 * t243 + t226 * t233;
t228 = sin(qJ(5));
t231 = cos(qJ(5));
t237 = qJD(2) * (t232 * t228 - t229 * t231);
t221 = t223 * qJD(2);
t219 = t222 * qJD(2);
t1 = [0, 0, t221, 0, -t221 (-t224 * t242 - t226 * t230) * t237 - t250 * ((-t223 * t229 + t224 * t244) * t228 - (t223 * t232 + t224 * t246) * t231); 0, 0, t219, 0, -t219 (-t224 * t230 + t226 * t242) * t237 + t250 * ((t222 * t229 + t226 * t244) * t228 + (t222 * t232 - t226 * t246) * t231); 0, 0, t241, 0, -t241, t233 * t225 * t237 + t250 * ((-t227 * t232 + t229 * t245) * t228 + (t227 * t229 + t230 * t244) * t231);];
JgD_rot  = t1;
