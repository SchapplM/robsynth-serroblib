% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR3_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:24
% EndTime: 2019-02-26 20:05:24
% DurationCPUTime: 0.08s
% Computational Cost: add. (39->18), mult. (136->43), div. (0->0), fcn. (146->12), ass. (0->28)
t211 = cos(pkin(7));
t205 = sin(pkin(13));
t209 = cos(pkin(13));
t213 = sin(qJ(3));
t215 = cos(qJ(3));
t204 = -t215 * t205 - t213 * t209;
t217 = qJD(3) * t204;
t199 = t211 * t217;
t227 = qJD(2) * t204 + t199;
t207 = sin(pkin(7));
t198 = t207 * t217;
t208 = sin(pkin(6));
t226 = t208 * t198;
t212 = cos(pkin(6));
t214 = sin(qJ(2));
t225 = t212 * t214;
t216 = cos(qJ(2));
t224 = t212 * t216;
t222 = t205 * t213 - t209 * t215;
t206 = sin(pkin(12));
t210 = cos(pkin(12));
t220 = t206 * t216 + t210 * t225;
t218 = t206 * t225 - t210 * t216;
t203 = t222 * qJD(3);
t202 = t218 * qJD(2);
t201 = t220 * qJD(2);
t200 = t222 * t211;
t1 = [0, 0, -t202 * t207, 0, t202 * t200 + t218 * t203 - t206 * t226 + t227 * (t206 * t224 + t210 * t214) 0; 0, 0, t201 * t207, 0, -t201 * t200 - t220 * t203 + t210 * t226 - t227 * (-t206 * t214 + t210 * t224) 0; 0, 0, t208 * qJD(2) * t214 * t207, 0, -t212 * t198 + (-t199 * t216 - t203 * t214 + (-t200 * t214 - t204 * t216) * qJD(2)) * t208, 0;];
JgD_rot  = t1;
