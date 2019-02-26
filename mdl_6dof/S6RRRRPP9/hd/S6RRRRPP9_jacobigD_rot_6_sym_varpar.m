% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPP9_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:16
% EndTime: 2019-02-26 22:30:16
% DurationCPUTime: 0.05s
% Computational Cost: add. (22->16), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->21)
t182 = sin(pkin(6));
t187 = cos(qJ(3));
t201 = t182 * t187;
t189 = cos(qJ(1));
t200 = t182 * t189;
t185 = sin(qJ(2));
t186 = sin(qJ(1));
t199 = t185 * t186;
t198 = t185 * t189;
t188 = cos(qJ(2));
t197 = t186 * t188;
t196 = t189 * t188;
t195 = qJD(1) * t182;
t184 = sin(qJ(3));
t194 = qJD(2) * t184;
t183 = cos(pkin(6));
t193 = t183 * t196 - t199;
t192 = t183 * t197 + t198;
t191 = t183 * t198 + t197;
t190 = -t183 * t199 + t196;
t1 = [0, t189 * t195, t193 * qJD(1) + t190 * qJD(2) (t186 * t182 * t184 + t190 * t187) * qJD(3) - t192 * t194 + (-t191 * t184 - t187 * t200) * qJD(1), 0, 0; 0, t186 * t195, t192 * qJD(1) + t191 * qJD(2) (-t184 * t200 + t191 * t187) * qJD(3) + t193 * t194 + (t190 * t184 - t186 * t201) * qJD(1), 0, 0; 0, 0, t182 * qJD(2) * t185, t182 * t188 * t194 + (t183 * t184 + t185 * t201) * qJD(3), 0, 0;];
JgD_rot  = t1;
