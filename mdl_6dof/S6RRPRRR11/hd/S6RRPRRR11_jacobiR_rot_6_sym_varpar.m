% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:46
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR11_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:45:55
% EndTime: 2019-02-22 11:45:55
% DurationCPUTime: 0.04s
% Computational Cost: add. (110->15), mult. (68->20), div. (0->0), fcn. (117->6), ass. (0->17)
t90 = sin(qJ(2));
t91 = sin(qJ(1));
t97 = t91 * t90;
t89 = qJ(4) + qJ(5) + qJ(6);
t87 = sin(t89);
t92 = cos(qJ(2));
t86 = t92 * t87;
t88 = cos(t89);
t96 = t92 * t88;
t93 = cos(qJ(1));
t95 = t93 * t90;
t94 = t93 * t92;
t85 = -t87 * t97 + t93 * t88;
t84 = t93 * t87 + t88 * t97;
t83 = t87 * t95 + t91 * t88;
t82 = -t91 * t87 + t88 * t95;
t1 = [t85, t87 * t94, 0, t82, t82, t82; t83, t91 * t86, 0, t84, t84, t84; 0, t90 * t87, 0, -t96, -t96, -t96; -t84, t88 * t94, 0, -t83, -t83, -t83; t82, t91 * t96, 0, t85, t85, t85; 0, t90 * t88, 0, t86, t86, t86; -t91 * t92, -t95, 0, 0, 0, 0; t94, -t97, 0, 0, 0, 0; 0, t92, 0, 0, 0, 0;];
JR_rot  = t1;
