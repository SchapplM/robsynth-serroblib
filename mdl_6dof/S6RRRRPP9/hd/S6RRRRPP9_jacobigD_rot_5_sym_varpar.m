% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JgD_rot = S6RRRRPP9_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobigD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:06
% EndTime: 2019-02-26 22:30:06
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->16), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->21)
t185 = sin(pkin(6));
t190 = cos(qJ(3));
t204 = t185 * t190;
t192 = cos(qJ(1));
t203 = t185 * t192;
t188 = sin(qJ(2));
t189 = sin(qJ(1));
t202 = t188 * t189;
t201 = t188 * t192;
t191 = cos(qJ(2));
t200 = t189 * t191;
t199 = t192 * t191;
t198 = qJD(1) * t185;
t187 = sin(qJ(3));
t197 = qJD(2) * t187;
t186 = cos(pkin(6));
t196 = t186 * t199 - t202;
t195 = t186 * t200 + t201;
t194 = t186 * t201 + t200;
t193 = -t186 * t202 + t199;
t1 = [0, t192 * t198, t196 * qJD(1) + t193 * qJD(2) (t189 * t185 * t187 + t193 * t190) * qJD(3) - t195 * t197 + (-t194 * t187 - t190 * t203) * qJD(1), 0, 0; 0, t189 * t198, t195 * qJD(1) + t194 * qJD(2) (-t187 * t203 + t194 * t190) * qJD(3) + t196 * t197 + (t193 * t187 - t189 * t204) * qJD(1), 0, 0; 0, 0, t185 * qJD(2) * t188, t185 * t191 * t197 + (t186 * t187 + t188 * t204) * qJD(3), 0, 0;];
JgD_rot  = t1;
