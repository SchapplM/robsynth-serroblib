% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JgD_rot = S6RRPRRP5_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:21
% EndTime: 2019-02-26 21:48:21
% DurationCPUTime: 0.09s
% Computational Cost: add. (45->22), mult. (147->54), div. (0->0), fcn. (162->10), ass. (0->24)
t204 = sin(pkin(11));
t206 = cos(pkin(11));
t209 = sin(qJ(2));
t212 = cos(qJ(2));
t215 = t209 * t204 - t212 * t206;
t222 = t215 * qJD(2);
t216 = t212 * t204 + t209 * t206;
t201 = t216 * qJD(2);
t205 = sin(pkin(6));
t210 = sin(qJ(1));
t221 = t205 * t210;
t213 = cos(qJ(1));
t220 = t205 * t213;
t219 = qJD(1) * t205;
t207 = cos(pkin(6));
t199 = t216 * t207;
t218 = t213 * t199 - t210 * t215;
t217 = -t210 * t199 - t213 * t215;
t211 = cos(qJ(4));
t208 = sin(qJ(4));
t198 = t215 * t207;
t197 = t207 * t222;
t196 = t207 * t201;
t1 = [0, t213 * t219, 0, -t210 * t196 - t213 * t222 + (-t198 * t213 - t210 * t216) * qJD(1) (t210 * t197 - t213 * t201) * t208 + (t208 * t221 + t217 * t211) * qJD(4) + (-t218 * t208 - t211 * t220) * qJD(1), 0; 0, t210 * t219, 0, t213 * t196 - t210 * t222 + (-t198 * t210 + t213 * t216) * qJD(1) (-t213 * t197 - t210 * t201) * t208 + (-t208 * t220 + t218 * t211) * qJD(4) + (t217 * t208 - t211 * t221) * qJD(1), 0; 0, 0, 0, t205 * t201, t207 * qJD(4) * t208 + (t216 * qJD(4) * t211 - t208 * t222) * t205, 0;];
JgD_rot  = t1;
