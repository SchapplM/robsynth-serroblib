% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR14V3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobigD_rot_6_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:09
% EndTime: 2019-04-12 15:12:09
% DurationCPUTime: 0.11s
% Computational Cost: add. (33->25), mult. (110->53), div. (0->0), fcn. (112->8), ass. (0->25)
t181 = sin(qJ(4));
t186 = cos(qJ(2));
t204 = -qJD(4) * t186 + qJD(1);
t205 = t204 * t181;
t180 = sin(qJ(5));
t185 = cos(qJ(4));
t203 = t180 * t185;
t202 = t185 * t186;
t183 = sin(qJ(1));
t201 = qJD(1) * t183;
t187 = cos(qJ(1));
t200 = qJD(1) * t187;
t182 = sin(qJ(2));
t199 = qJD(2) * t182;
t198 = qJD(2) * t185;
t197 = qJD(2) * t186;
t196 = qJD(4) * t182;
t194 = t187 * qJD(2);
t192 = qJD(1) * t186 - qJD(4);
t191 = -qJD(5) + t198;
t190 = t204 * t185;
t189 = t192 * t183;
t188 = -t182 * t201 + t186 * t194;
t184 = cos(qJ(5));
t1 = [0, t200, 0, t188, -t187 * t190 + (-t182 * t194 - t189) * t181, -t189 * t203 + (-t191 * t182 + t205) * t180 * t187 + ((t183 * t181 + t187 * t202) * qJD(5) - t188) * t184; 0, t201, 0, t182 * t200 + t183 * t197, -t183 * t190 + (-t183 * t199 + t192 * t187) * t181 (t192 * t203 + (-qJD(1) * t182 - t181 * qJD(5)) * t184) * t187 + ((-t182 * t198 + t205) * t180 - t184 * t197 + (t182 * t180 + t184 * t202) * qJD(5)) * t183; 0, 0, 0, t199, t181 * t197 + t185 * t196 (qJD(5) * t185 - qJD(2)) * t184 * t182 + (-t181 * t196 + t191 * t186) * t180;];
JgD_rot  = t1;
