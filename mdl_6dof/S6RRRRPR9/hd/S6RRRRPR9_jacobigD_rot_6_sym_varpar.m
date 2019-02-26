% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR9_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:13
% EndTime: 2019-02-26 22:35:13
% DurationCPUTime: 0.07s
% Computational Cost: add. (45->18), mult. (100->40), div. (0->0), fcn. (102->8), ass. (0->27)
t203 = sin(pkin(6));
t206 = sin(qJ(1));
t221 = t203 * t206;
t208 = cos(qJ(1));
t220 = t203 * t208;
t205 = sin(qJ(2));
t219 = t205 * t206;
t218 = t205 * t208;
t207 = cos(qJ(2));
t217 = t206 * t207;
t216 = t208 * t207;
t215 = qJD(1) * t203;
t202 = qJ(3) + qJ(4);
t199 = sin(t202);
t214 = qJD(2) * t199;
t213 = qJD(2) * t203;
t204 = cos(pkin(6));
t212 = t204 * t216 - t219;
t211 = t204 * t217 + t218;
t210 = t204 * t218 + t217;
t209 = -t204 * t219 + t216;
t201 = qJD(3) + qJD(4);
t200 = cos(t202);
t198 = t205 * t213;
t197 = t211 * qJD(1) + t210 * qJD(2);
t196 = t212 * qJD(1) + t209 * qJD(2);
t1 = [0, t208 * t215, t196, t196, 0 (t199 * t221 + t209 * t200) * t201 - t211 * t214 + (-t210 * t199 - t200 * t220) * qJD(1); 0, t206 * t215, t197, t197, 0 (-t199 * t220 + t210 * t200) * t201 + t212 * t214 + (t209 * t199 - t200 * t221) * qJD(1); 0, 0, t198, t198, 0, t203 * t205 * t201 * t200 + (t201 * t204 + t207 * t213) * t199;];
JgD_rot  = t1;
