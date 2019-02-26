% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRRR10_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:50
% EndTime: 2019-02-26 21:19:50
% DurationCPUTime: 0.07s
% Computational Cost: add. (44->18), mult. (166->42), div. (0->0), fcn. (174->10), ass. (0->29)
t210 = sin(pkin(7));
t211 = sin(pkin(6));
t231 = t211 * t210;
t213 = cos(pkin(7));
t217 = cos(qJ(3));
t230 = t213 * t217;
t209 = sin(pkin(13));
t216 = sin(qJ(1));
t229 = t216 * t209;
t212 = cos(pkin(13));
t228 = t216 * t212;
t218 = cos(qJ(1));
t227 = t218 * t209;
t226 = t218 * t212;
t225 = t216 * t231;
t224 = t218 * t231;
t223 = qJD(1) * t211 * t213;
t214 = cos(pkin(6));
t222 = t214 * t226 - t229;
t221 = -t214 * t228 - t227;
t220 = t214 * t227 + t228;
t219 = -t214 * t229 + t226;
t215 = sin(qJ(3));
t208 = t221 * qJD(1);
t207 = t222 * qJD(1);
t206 = (t210 * t214 * t215 + (t212 * t213 * t215 + t209 * t217) * t211) * qJD(3);
t205 = -t208 * t230 + (t220 * t217 + (t222 * t213 - t224) * t215) * qJD(3) + (t219 * t215 - t217 * t225) * qJD(1);
t204 = t207 * t230 + (t219 * t217 + (t221 * t213 + t225) * t215) * qJD(3) + (-t220 * t215 - t217 * t224) * qJD(1);
t1 = [0, 0, t207 * t210 + t218 * t223, t204, t204, 0; 0, 0, -t208 * t210 + t216 * t223, t205, t205, 0; 0, 0, 0, t206, t206, 0;];
JgD_rot  = t1;
