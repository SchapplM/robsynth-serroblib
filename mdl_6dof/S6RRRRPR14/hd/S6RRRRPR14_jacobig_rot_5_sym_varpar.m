% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPR14_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:12
% EndTime: 2019-02-26 22:38:12
% DurationCPUTime: 0.04s
% Computational Cost: add. (16->14), mult. (48->31), div. (0->0), fcn. (72->10), ass. (0->19)
t179 = sin(pkin(6));
t184 = sin(qJ(1));
t193 = t184 * t179;
t183 = sin(qJ(2));
t192 = t184 * t183;
t186 = cos(qJ(2));
t191 = t184 * t186;
t187 = cos(qJ(1));
t190 = t187 * t179;
t189 = t187 * t183;
t188 = t187 * t186;
t185 = cos(qJ(3));
t182 = sin(qJ(3));
t181 = cos(pkin(6));
t180 = cos(pkin(7));
t178 = sin(pkin(7));
t177 = -t181 * t191 - t189;
t176 = t181 * t188 - t192;
t1 = [0, t193, -t177 * t178 + t180 * t193 (-t181 * t192 + t188) * t182 + (-t177 * t180 - t178 * t193) * t185, 0, 0; 0, -t190, -t176 * t178 - t180 * t190 (t181 * t189 + t191) * t182 + (-t176 * t180 + t178 * t190) * t185, 0, 0; 1, t181, -t179 * t186 * t178 + t181 * t180, -t181 * t178 * t185 + (-t180 * t185 * t186 + t182 * t183) * t179, 0, 0;];
Jg_rot  = t1;
