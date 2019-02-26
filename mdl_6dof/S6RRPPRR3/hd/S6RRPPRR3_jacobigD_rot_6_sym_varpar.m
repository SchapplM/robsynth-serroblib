% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPPRR3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:34
% EndTime: 2019-02-26 21:29:34
% DurationCPUTime: 0.09s
% Computational Cost: add. (56->23), mult. (147->54), div. (0->0), fcn. (162->10), ass. (0->25)
t212 = sin(pkin(11));
t214 = cos(pkin(11));
t216 = sin(qJ(2));
t218 = cos(qJ(2));
t221 = t216 * t212 - t218 * t214;
t228 = t221 * qJD(2);
t222 = t218 * t212 + t216 * t214;
t206 = t222 * qJD(2);
t213 = sin(pkin(6));
t217 = sin(qJ(1));
t227 = t213 * t217;
t219 = cos(qJ(1));
t226 = t213 * t219;
t225 = qJD(1) * t213;
t215 = cos(pkin(6));
t204 = t222 * t215;
t224 = t219 * t204 - t217 * t221;
t223 = -t217 * t204 - t219 * t221;
t211 = pkin(12) + qJ(5);
t210 = cos(t211);
t209 = sin(t211);
t203 = t221 * t215;
t202 = t215 * t228;
t201 = t215 * t206;
t1 = [0, t219 * t225, 0, 0, -t217 * t201 - t219 * t228 + (-t203 * t219 - t217 * t222) * qJD(1) (t217 * t202 - t219 * t206) * t209 + (t209 * t227 + t223 * t210) * qJD(5) + (-t224 * t209 - t210 * t226) * qJD(1); 0, t217 * t225, 0, 0, t219 * t201 - t217 * t228 + (-t203 * t217 + t219 * t222) * qJD(1) (-t219 * t202 - t217 * t206) * t209 + (-t209 * t226 + t224 * t210) * qJD(5) + (t223 * t209 - t210 * t227) * qJD(1); 0, 0, 0, 0, t213 * t206, t215 * qJD(5) * t209 + (t222 * qJD(5) * t210 - t209 * t228) * t213;];
JgD_rot  = t1;
