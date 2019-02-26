% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR6
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
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR6_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:05
% EndTime: 2019-02-26 20:07:05
% DurationCPUTime: 0.06s
% Computational Cost: add. (24->17), mult. (88->43), div. (0->0), fcn. (92->10), ass. (0->24)
t190 = sin(pkin(7));
t191 = sin(pkin(6));
t209 = t191 * t190;
t193 = cos(pkin(7));
t197 = cos(qJ(3));
t208 = t193 * t197;
t194 = cos(pkin(6));
t196 = sin(qJ(2));
t207 = t194 * t196;
t198 = cos(qJ(2));
t206 = t194 * t198;
t195 = sin(qJ(3));
t205 = t195 * t198;
t204 = t196 * t197;
t203 = qJD(2) * t195;
t189 = sin(pkin(12));
t192 = cos(pkin(12));
t202 = -t189 * t196 + t192 * t206;
t201 = t189 * t198 + t192 * t207;
t200 = -t189 * t206 - t192 * t196;
t199 = t189 * t207 - t192 * t198;
t188 = t199 * qJD(2);
t187 = t201 * qJD(2);
t1 = [0, 0, -t188 * t190, 0, -t188 * t208 + t200 * t203 + (-t199 * t197 + (t189 * t209 + t200 * t193) * t195) * qJD(3), 0; 0, 0, t187 * t190, 0, t187 * t208 + t202 * t203 + (t201 * t197 + (-t192 * t209 + t202 * t193) * t195) * qJD(3), 0; 0, 0, qJD(2) * t196 * t209, 0, t194 * t190 * qJD(3) * t195 + ((t193 * t205 + t204) * qJD(3) + (t193 * t204 + t205) * qJD(2)) * t191, 0;];
JgD_rot  = t1;
