% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRPR12_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobigD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:17
% EndTime: 2019-02-26 21:07:17
% DurationCPUTime: 0.07s
% Computational Cost: add. (24->18), mult. (92->42), div. (0->0), fcn. (96->10), ass. (0->26)
t205 = sin(pkin(7));
t206 = sin(pkin(6));
t226 = t206 * t205;
t208 = cos(pkin(7));
t212 = cos(qJ(3));
t225 = t208 * t212;
t204 = sin(pkin(12));
t211 = sin(qJ(1));
t224 = t211 * t204;
t207 = cos(pkin(12));
t223 = t211 * t207;
t213 = cos(qJ(1));
t222 = t213 * t204;
t221 = t213 * t207;
t220 = t211 * t226;
t219 = t213 * t226;
t218 = qJD(1) * t206 * t208;
t209 = cos(pkin(6));
t217 = t209 * t221 - t224;
t216 = -t209 * t223 - t222;
t215 = t209 * t222 + t223;
t214 = -t209 * t224 + t221;
t210 = sin(qJ(3));
t203 = t216 * qJD(1);
t202 = t217 * qJD(1);
t1 = [0, 0, t202 * t205 + t213 * t218, t202 * t225 + (t214 * t212 + (t216 * t208 + t220) * t210) * qJD(3) + (-t215 * t210 - t212 * t219) * qJD(1), 0, 0; 0, 0, -t203 * t205 + t211 * t218, -t203 * t225 + (t215 * t212 + (t217 * t208 - t219) * t210) * qJD(3) + (t214 * t210 - t212 * t220) * qJD(1), 0, 0; 0, 0, 0 (t205 * t209 * t210 + (t207 * t208 * t210 + t204 * t212) * t206) * qJD(3), 0, 0;];
JgD_rot  = t1;
