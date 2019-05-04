% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR10V2_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobig_rot_6_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:32
% EndTime: 2019-04-11 14:56:32
% DurationCPUTime: 0.03s
% Computational Cost: add. (22->11), mult. (24->18), div. (0->0), fcn. (44->8), ass. (0->16)
t154 = qJ(2) + qJ(3);
t152 = sin(t154);
t157 = sin(qJ(1));
t166 = t157 * t152;
t156 = sin(qJ(4));
t165 = t157 * t156;
t159 = cos(qJ(4));
t164 = t157 * t159;
t160 = cos(qJ(1));
t163 = t160 * t152;
t162 = t160 * t156;
t161 = t160 * t159;
t158 = cos(qJ(5));
t155 = sin(qJ(5));
t153 = cos(t154);
t1 = [0, t157, t157, t163, t153 * t162 - t164 (t153 * t161 + t165) * t155 - t158 * t163; 0, -t160, -t160, t166, t153 * t165 + t161 (t153 * t164 - t162) * t155 - t158 * t166; 1, 0, 0, -t153, t152 * t156, t152 * t159 * t155 + t153 * t158;];
Jg_rot  = t1;
