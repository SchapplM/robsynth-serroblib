% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPPRRR1_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobigD_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:44
% EndTime: 2019-02-26 19:38:44
% DurationCPUTime: 0.08s
% Computational Cost: add. (31->23), mult. (101->55), div. (0->0), fcn. (126->14), ass. (0->27)
t208 = sin(pkin(12));
t217 = cos(pkin(6));
t227 = t208 * t217;
t210 = sin(pkin(7));
t211 = sin(pkin(6));
t226 = t210 * t211;
t225 = t210 * t217;
t216 = cos(pkin(7));
t224 = t211 * t216;
t213 = cos(pkin(13));
t223 = t213 * t216;
t214 = cos(pkin(12));
t222 = t214 * t217;
t207 = sin(pkin(13));
t202 = -t208 * t207 + t213 * t222;
t221 = t202 * t216 - t214 * t226;
t204 = -t214 * t207 - t213 * t227;
t220 = t204 * t216 + t208 * t226;
t219 = cos(qJ(4));
t218 = sin(qJ(4));
t215 = cos(pkin(8));
t212 = cos(pkin(14));
t209 = sin(pkin(8));
t206 = sin(pkin(14));
t205 = -t207 * t227 + t214 * t213;
t203 = t207 * t222 + t208 * t213;
t1 = [0, 0, 0, 0 ((t205 * t212 + t220 * t206) * t219 + ((-t205 * t206 + t220 * t212) * t215 + (-t204 * t210 + t208 * t224) * t209) * t218) * qJD(4), 0; 0, 0, 0, 0 ((t203 * t212 + t221 * t206) * t219 + ((-t203 * t206 + t221 * t212) * t215 + (-t202 * t210 - t214 * t224) * t209) * t218) * qJD(4), 0; 0, 0, 0, 0 ((t211 * t207 * t212 + (t211 * t223 + t225) * t206) * t219 + ((t212 * t225 + (-t206 * t207 + t212 * t223) * t211) * t215 + (-t213 * t226 + t217 * t216) * t209) * t218) * qJD(4), 0;];
JgD_rot  = t1;
