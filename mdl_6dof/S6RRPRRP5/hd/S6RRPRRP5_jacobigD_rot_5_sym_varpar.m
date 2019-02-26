% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRP5_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_jacobigD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:21
% EndTime: 2019-02-26 21:48:21
% DurationCPUTime: 0.09s
% Computational Cost: add. (45->22), mult. (147->54), div. (0->0), fcn. (162->10), ass. (0->24)
t202 = sin(pkin(11));
t204 = cos(pkin(11));
t207 = sin(qJ(2));
t210 = cos(qJ(2));
t213 = t207 * t202 - t210 * t204;
t220 = t213 * qJD(2);
t214 = t210 * t202 + t207 * t204;
t199 = t214 * qJD(2);
t203 = sin(pkin(6));
t208 = sin(qJ(1));
t219 = t203 * t208;
t211 = cos(qJ(1));
t218 = t203 * t211;
t217 = qJD(1) * t203;
t205 = cos(pkin(6));
t197 = t214 * t205;
t216 = t211 * t197 - t208 * t213;
t215 = -t208 * t197 - t211 * t213;
t209 = cos(qJ(4));
t206 = sin(qJ(4));
t196 = t213 * t205;
t195 = t205 * t220;
t194 = t205 * t199;
t1 = [0, t211 * t217, 0, -t208 * t194 - t211 * t220 + (-t196 * t211 - t208 * t214) * qJD(1) (t208 * t195 - t211 * t199) * t206 + (t206 * t219 + t215 * t209) * qJD(4) + (-t216 * t206 - t209 * t218) * qJD(1), 0; 0, t208 * t217, 0, t211 * t194 - t208 * t220 + (-t196 * t208 + t211 * t214) * qJD(1) (-t211 * t195 - t208 * t199) * t206 + (-t206 * t218 + t216 * t209) * qJD(4) + (t215 * t206 - t209 * t219) * qJD(1), 0; 0, 0, 0, t203 * t199, t205 * qJD(4) * t206 + (t214 * qJD(4) * t209 - t206 * t220) * t203, 0;];
JgD_rot  = t1;
