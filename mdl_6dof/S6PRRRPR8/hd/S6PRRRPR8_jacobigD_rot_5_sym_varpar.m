% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR8_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobigD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:42
% EndTime: 2019-02-26 20:14:42
% DurationCPUTime: 0.08s
% Computational Cost: add. (24->17), mult. (88->43), div. (0->0), fcn. (92->10), ass. (0->24)
t203 = sin(pkin(7));
t204 = sin(pkin(6));
t222 = t204 * t203;
t206 = cos(pkin(7));
t210 = cos(qJ(3));
t221 = t206 * t210;
t207 = cos(pkin(6));
t209 = sin(qJ(2));
t220 = t207 * t209;
t211 = cos(qJ(2));
t219 = t207 * t211;
t208 = sin(qJ(3));
t218 = t208 * t211;
t217 = t209 * t210;
t216 = qJD(2) * t208;
t202 = sin(pkin(12));
t205 = cos(pkin(12));
t215 = -t202 * t209 + t205 * t219;
t214 = t202 * t211 + t205 * t220;
t213 = -t202 * t219 - t205 * t209;
t212 = t202 * t220 - t205 * t211;
t201 = t212 * qJD(2);
t200 = t214 * qJD(2);
t1 = [0, 0, -t201 * t203, -t201 * t221 + t213 * t216 + (-t212 * t210 + (t202 * t222 + t213 * t206) * t208) * qJD(3), 0, 0; 0, 0, t200 * t203, t200 * t221 + t215 * t216 + (t214 * t210 + (-t205 * t222 + t215 * t206) * t208) * qJD(3), 0, 0; 0, 0, qJD(2) * t209 * t222, t207 * t203 * qJD(3) * t208 + ((t206 * t218 + t217) * qJD(3) + (t206 * t217 + t218) * qJD(2)) * t204, 0, 0;];
JgD_rot  = t1;
