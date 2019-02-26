% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRPR4_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:39:33
% EndTime: 2019-02-26 21:39:33
% DurationCPUTime: 0.09s
% Computational Cost: add. (56->23), mult. (147->54), div. (0->0), fcn. (162->10), ass. (0->25)
t213 = sin(pkin(11));
t215 = cos(pkin(11));
t217 = sin(qJ(2));
t219 = cos(qJ(2));
t222 = t217 * t213 - t219 * t215;
t229 = t222 * qJD(2);
t223 = t219 * t213 + t217 * t215;
t207 = t223 * qJD(2);
t214 = sin(pkin(6));
t218 = sin(qJ(1));
t228 = t214 * t218;
t220 = cos(qJ(1));
t227 = t214 * t220;
t226 = qJD(1) * t214;
t216 = cos(pkin(6));
t205 = t223 * t216;
t225 = t220 * t205 - t218 * t222;
t224 = -t218 * t205 - t220 * t222;
t212 = qJ(4) + pkin(12);
t211 = cos(t212);
t210 = sin(t212);
t204 = t222 * t216;
t203 = t216 * t229;
t202 = t216 * t207;
t1 = [0, t220 * t226, 0, -t218 * t202 - t220 * t229 + (-t204 * t220 - t218 * t223) * qJD(1), 0 (t218 * t203 - t220 * t207) * t210 + (t210 * t228 + t224 * t211) * qJD(4) + (-t225 * t210 - t211 * t227) * qJD(1); 0, t218 * t226, 0, t220 * t202 - t218 * t229 + (-t204 * t218 + t220 * t223) * qJD(1), 0 (-t220 * t203 - t218 * t207) * t210 + (-t210 * t227 + t225 * t211) * qJD(4) + (t224 * t210 - t211 * t228) * qJD(1); 0, 0, 0, t214 * t207, 0, t216 * qJD(4) * t210 + (t223 * qJD(4) * t211 - t210 * t229) * t214;];
JgD_rot  = t1;
