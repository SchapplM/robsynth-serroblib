% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPPRR4_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:11
% EndTime: 2019-02-26 21:30:11
% DurationCPUTime: 0.09s
% Computational Cost: add. (45->23), mult. (147->54), div. (0->0), fcn. (162->10), ass. (0->24)
t203 = sin(pkin(11));
t205 = cos(pkin(11));
t208 = sin(qJ(2));
t211 = cos(qJ(2));
t214 = t208 * t203 - t211 * t205;
t221 = t214 * qJD(2);
t215 = t211 * t203 + t208 * t205;
t200 = t215 * qJD(2);
t204 = sin(pkin(6));
t209 = sin(qJ(1));
t220 = t204 * t209;
t212 = cos(qJ(1));
t219 = t204 * t212;
t218 = qJD(1) * t204;
t206 = cos(pkin(6));
t197 = t214 * t206;
t217 = -t212 * t197 - t209 * t215;
t216 = -t209 * t197 + t212 * t215;
t210 = cos(qJ(5));
t207 = sin(qJ(5));
t198 = t215 * t206;
t196 = t206 * t221;
t195 = t206 * t200;
t1 = [0, t212 * t218, 0, 0, t209 * t196 - t212 * t200 + (-t198 * t212 + t209 * t214) * qJD(1) -(-t209 * t195 - t212 * t221) * t210 + (t216 * t207 + t210 * t220) * qJD(5) + (t207 * t219 - t217 * t210) * qJD(1); 0, t209 * t218, 0, 0, -t212 * t196 - t209 * t200 + (-t198 * t209 - t212 * t214) * qJD(1) -(t212 * t195 - t209 * t221) * t210 + (-t217 * t207 - t210 * t219) * qJD(5) + (t207 * t220 - t216 * t210) * qJD(1); 0, 0, 0, 0, -t204 * t221, t206 * qJD(5) * t210 + (t214 * qJD(5) * t207 - t210 * t200) * t204;];
JgD_rot  = t1;
