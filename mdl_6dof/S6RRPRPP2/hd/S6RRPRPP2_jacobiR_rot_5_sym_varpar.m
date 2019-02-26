% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPP2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:21
% EndTime: 2019-02-26 21:35:21
% DurationCPUTime: 0.03s
% Computational Cost: add. (35->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
t88 = sin(qJ(4));
t89 = sin(qJ(1));
t95 = t89 * t88;
t90 = cos(qJ(4));
t94 = t89 * t90;
t91 = cos(qJ(1));
t93 = t91 * t88;
t92 = t91 * t90;
t87 = qJ(2) + pkin(9);
t86 = cos(t87);
t85 = sin(t87);
t84 = t86 * t92 + t95;
t83 = t86 * t93 - t94;
t82 = t86 * t94 - t93;
t81 = -t86 * t95 - t92;
t1 = [-t82, -t85 * t92, 0, -t83, 0, 0; t84, -t85 * t94, 0, t81, 0, 0; 0, t86 * t90, 0, -t85 * t88, 0, 0; -t89 * t85, t91 * t86, 0, 0, 0, 0; t91 * t85, t89 * t86, 0, 0, 0, 0; 0, t85, 0, 0, 0, 0; t81, -t85 * t93, 0, t84, 0, 0; t83, -t85 * t95, 0, t82, 0, 0; 0, t86 * t88, 0, t85 * t90, 0, 0;];
JR_rot  = t1;
