% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRPRR9_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:24
% EndTime: 2019-02-26 20:53:24
% DurationCPUTime: 0.07s
% Computational Cost: add. (39->20), mult. (141->45), div. (0->0), fcn. (151->12), ass. (0->34)
t212 = sin(pkin(6));
t218 = sin(qJ(1));
t233 = t212 * t218;
t220 = cos(qJ(1));
t232 = t212 * t220;
t210 = sin(pkin(12));
t231 = t218 * t210;
t214 = cos(pkin(12));
t230 = t218 * t214;
t229 = t220 * t210;
t228 = t220 * t214;
t215 = cos(pkin(7));
t227 = qJD(1) * t212 * t215;
t209 = sin(pkin(13));
t213 = cos(pkin(13));
t217 = sin(qJ(3));
t219 = cos(qJ(3));
t208 = -t219 * t209 - t217 * t213;
t226 = t209 * t217 - t213 * t219;
t216 = cos(pkin(6));
t225 = t216 * t228 - t231;
t224 = -t216 * t230 - t229;
t223 = t216 * t229 + t230;
t222 = t216 * t231 - t228;
t221 = qJD(3) * t208;
t211 = sin(pkin(7));
t207 = t226 * qJD(3);
t206 = t224 * qJD(1);
t205 = t225 * qJD(1);
t204 = t226 * t215;
t203 = t226 * t211;
t202 = t215 * t221;
t201 = t211 * t221;
t1 = [0, 0, t205 * t211 + t220 * t227, 0, t222 * t207 - t205 * t204 - t224 * t202 - t201 * t233 + (t203 * t232 + t223 * t208) * qJD(1), 0; 0, 0, -t206 * t211 + t218 * t227, 0, -t223 * t207 + t206 * t204 - t225 * t202 + t201 * t232 + (t203 * t233 + t222 * t208) * qJD(1), 0; 0, 0, 0, 0, -t216 * t201 + (-t202 * t214 - t207 * t210) * t212, 0;];
JgD_rot  = t1;
