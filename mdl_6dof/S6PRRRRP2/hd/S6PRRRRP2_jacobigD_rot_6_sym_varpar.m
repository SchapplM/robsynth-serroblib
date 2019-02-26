% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRP2_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:43
% EndTime: 2019-02-26 20:15:43
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->12), mult. (60->30), div. (0->0), fcn. (62->8), ass. (0->21)
t210 = qJ(3) + qJ(4);
t207 = sin(t210);
t212 = sin(pkin(6));
t223 = t212 * t207;
t214 = cos(pkin(6));
t215 = sin(qJ(2));
t222 = t214 * t215;
t216 = cos(qJ(2));
t221 = t214 * t216;
t220 = qJD(2) * t207;
t219 = qJD(2) * t212;
t211 = sin(pkin(11));
t213 = cos(pkin(11));
t218 = t211 * t216 + t213 * t222;
t217 = -t211 * t222 + t213 * t216;
t209 = qJD(3) + qJD(4);
t208 = cos(t210);
t206 = t215 * t219;
t205 = t217 * qJD(2);
t204 = t218 * qJD(2);
t1 = [0, 0, t205, t205 (t217 * t208 + t211 * t223) * t209 + (-t211 * t221 - t213 * t215) * t220, 0; 0, 0, t204, t204 (t218 * t208 - t213 * t223) * t209 + (-t211 * t215 + t213 * t221) * t220, 0; 0, 0, t206, t206, t212 * t215 * t209 * t208 + (t209 * t214 + t216 * t219) * t207, 0;];
JgD_rot  = t1;
