% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR8_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:58
% EndTime: 2019-02-26 22:19:58
% DurationCPUTime: 0.06s
% Computational Cost: add. (60->17), mult. (130->40), div. (0->0), fcn. (134->8), ass. (0->26)
t190 = sin(pkin(6));
t193 = sin(qJ(1));
t208 = t190 * t193;
t195 = cos(qJ(1));
t207 = t190 * t195;
t192 = sin(qJ(2));
t206 = t192 * t193;
t205 = t192 * t195;
t194 = cos(qJ(2));
t204 = t193 * t194;
t203 = t195 * t194;
t202 = qJD(1) * t190;
t189 = qJ(3) + pkin(12);
t187 = sin(t189);
t201 = qJD(2) * t187;
t200 = qJD(2) * t190;
t191 = cos(pkin(6));
t199 = t191 * t203 - t206;
t198 = t191 * t204 + t205;
t197 = t191 * t205 + t204;
t196 = -t191 * t206 + t203;
t188 = cos(t189);
t186 = t194 * t187 * t200 + (t188 * t190 * t192 + t187 * t191) * qJD(3);
t185 = (-t187 * t207 + t197 * t188) * qJD(3) + t199 * t201 + (t196 * t187 - t188 * t208) * qJD(1);
t184 = (t187 * t208 + t196 * t188) * qJD(3) - t198 * t201 + (-t197 * t187 - t188 * t207) * qJD(1);
t1 = [0, t195 * t202, t199 * qJD(1) + t196 * qJD(2), 0, t184, t184; 0, t193 * t202, t198 * qJD(1) + t197 * qJD(2), 0, t185, t185; 0, 0, t192 * t200, 0, t186, t186;];
JgD_rot  = t1;
