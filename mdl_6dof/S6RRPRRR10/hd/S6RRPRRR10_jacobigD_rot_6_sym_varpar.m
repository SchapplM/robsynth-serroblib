% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR10_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:18
% EndTime: 2019-02-26 21:59:18
% DurationCPUTime: 0.06s
% Computational Cost: add. (60->17), mult. (130->40), div. (0->0), fcn. (134->8), ass. (0->26)
t189 = sin(pkin(6));
t192 = sin(qJ(1));
t207 = t189 * t192;
t194 = cos(qJ(1));
t206 = t189 * t194;
t191 = sin(qJ(2));
t205 = t191 * t192;
t204 = t191 * t194;
t193 = cos(qJ(2));
t203 = t192 * t193;
t202 = t194 * t193;
t201 = qJD(1) * t189;
t188 = pkin(12) + qJ(4);
t186 = sin(t188);
t200 = qJD(2) * t186;
t199 = qJD(2) * t189;
t190 = cos(pkin(6));
t198 = t190 * t202 - t205;
t197 = t190 * t203 + t204;
t196 = t190 * t204 + t203;
t195 = -t190 * t205 + t202;
t187 = cos(t188);
t185 = t193 * t186 * t199 + (t187 * t189 * t191 + t186 * t190) * qJD(4);
t184 = (-t186 * t206 + t196 * t187) * qJD(4) + t198 * t200 + (t195 * t186 - t187 * t207) * qJD(1);
t183 = (t186 * t207 + t195 * t187) * qJD(4) - t197 * t200 + (-t196 * t186 - t187 * t206) * qJD(1);
t1 = [0, t194 * t201, 0, t198 * qJD(1) + t195 * qJD(2), t183, t183; 0, t192 * t201, 0, t197 * qJD(1) + t196 * qJD(2), t184, t184; 0, 0, 0, t191 * t199, t185, t185;];
JgD_rot  = t1;
